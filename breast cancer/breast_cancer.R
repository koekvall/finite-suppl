set.seed(39)
library(discSurv)
library(penalized)
library(fsnet)
library(caret)
library(glmnet)
data(nki70, package = "penalized")

# Standardize non-clinical X for simplicity
nki70[, 8:77] <- scale(nki70[, 8:77])


###############################################################################
# Analysis with clinical predictors only
###############################################################################
library(splines2)
yl <- 3 * floor(nki70$time / 3) # Lower endpoint of 3-month interval
yu <- yl + 3 # Upper endpoint of interval
yu[nki70$event == 0] <- Inf # If censored, may have occurred later
times <- seq(1e-5, 72, by = 1e-5)
spline_mat <- iSpline(log(times), knots = NULL, degree = 3, intercept = FALSE)

X <- model.matrix(~., data = nki70[, 3:7])
Z <- cbind(matrix(0, 2 * nrow(X), ncol(spline_mat)), kronecker(X, -c(1, 1)))
for(ii in 1:nrow(X)){
  if((yl[ii] != 0)){
    Z[2 * (ii - 1) + 1, 1:ncol(spline_mat)] <- spline_mat[1e5 * yl[ii], ]
  } else{
    Z[2 * (ii - 1) + 1, ] <- 0
  }
  if(is.finite(yu[ii])){
    Z[2 * (ii - 1) + 2, 1:ncol(spline_mat)] <- spline_mat[1e5 * yu[ii], ]
  } else{
    Z[2 * (ii - 1) + 2, ] <- 0
  }
}

M <- log(cbind(yl, yu))
M[is.finite(M)] <- 0
box_constr <- cbind(rep(-Inf, ncol(Z)), rep(Inf, ncol(Z)))
box_constr[1:ncol(spline_mat), 1] <- 0 # Positive spline coefficients

flex_fit <- fsnet_flex(M = M, Z = Z,
                       theta = c(rep(0.1, 5), rep(0, ncol(Z) - 5)),
                       box_constr = box_constr,
                       verbose = FALSE,
                       method = "prox_newt",
                       distr = "ee",
                       lam = 0,
                       maxit = rep(1e4, 3))

simple_fit <- fsnet(Y = log(cbind(yl, yu)),
                    X = model.matrix(~., data = nki70[, 3:7]),
                    lam = 0,
                    distr = "ee",
                    fix_var = FALSE,
                    maxit = rep(1e4, 3),
                    tol = rep(1e-7, 2))

exp_fit <- fsnet(Y = log(cbind(yl, yu)),
                 X = model.matrix(~., data = nki70[, 3:7]),
                 lam = 0,
                 distr = "ee",
                 fix_var = TRUE,
                 maxit = rep(1e4, 3),
                 tol = rep(1e-8, 2))

null_fit <- fsnet(Y = log(cbind(yl, yu)),
                  X = matrix(1, ncol = 1, nrow = nrow(X)),
                  lam = 0,
                  distr = "ee",
                  fix_var = TRUE,
                  maxit = rep(1e4, 3),
                  tol = rep(1e-7, 2))

# Test models
pchisq(2 * (flex_fit$loglik - simple_fit$loglik), df = 2, lower = F)

pchisq(2 * (simple_fit$loglik - exp_fit$loglik), df = 1, lower = F)

pchisq(2 * (flex_fit$loglik - exp_fit$loglik), df = 1, lower = F)


# Standard errors and p-values

se <- sqrt(diag(solve(-hessian_fsnet(Y = log(cbind(yl, yu)),
                                     X = X,
                                     theta = exp_fit$theta[1:8],
                                     fix_s = TRUE,
                                     distr = "ee"))))
2 * pnorm(abs(exp_fit$theta[2:8] / se), lower = F)


# Estimated survival function at median predictors
Z_time <- cbind(spline_mat, -matrix(rep(apply(X, 2, median),
                                        each = nrow(spline_mat)),
                                    nrow = nrow(spline_mat)))

flex_surv <- exp(-exp(Z_time %*% flex_fit$theta))

Z_time_simp <- cbind(log(times), -matrix(rep(apply(X, 2, median),
                                             each = nrow(spline_mat)),
                                         nrow = nrow(spline_mat)))
theta_simp <- c(1, exp_fit$beta) /  exp_fit$sigma
exp_surv <- exp(-exp(Z_time_simp %*% theta_simp))

# Create figure

pdf("~/GitHub/finite-suppl/breast cancer/fig_bc.pdf", width = 12.5, height = 6)
par(cex.axis = 1.3, cex.lab = 1.3)
par(mfrow = c(1, 2))
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))

plot(x = times[floor(seq(1, length(times), length.out = 100))],
     y = flex_surv[floor(seq(1, length(times), length.out = 100))],
     type = "l",
     ylim = c(0, 1),
     lty = 1, lwd = 2,
     ylab = "Survival probability",
     xlab = "Months")

legend("topright", c("Flexible", "Exponential"), lty = 1:2, lwd = 2,
       bty = "n")

lines(x = times[floor(seq(1, length(times), length.out = 100))],
      y = exp_surv[floor(seq(1, length(times), length.out = 100))],
      lwd = 2,
      lty = 2)

barplot(rbind("Flexible" = flex_surv[1e5 * c(3, 6, 9, 12, 15)],
              "Exponential" = exp_surv[1e5 * c(3, 6, 9, 12, 15)]),
        beside = T, legend = T, names.arg = c("3", "6", "9", "12", "15"),
        ylim = c(0, 1),
        xlab = "Month",
        ylab = "Survival probability"
)

dev.off()



###############################################################################
# Analysis w. gene expressions
###############################################################################

# Split data into test and training set
train_idx <- rep(0, nrow(nki70))
train_idx[sample(1:nrow(nki70), 96)] <- 1
train_idx <- as.logical(train_idx)
train_dat <- nki70[train_idx, ]
test_dat <- nki70[!train_idx, ]

# Fold ID for glmnet cross-validation
train_dat$fold <- rep(1:10, length.out = nrow(train_dat))[sample(1:nrow(train_dat))]

# Fit our model
yl <- 3 * floor(train_dat$time / 3) # Lower endpoint of 3-month interval
yu <- yl + 3 # Upper endpoint of interval
yu[train_dat$event == 0] <- Inf # If censored, may have occurred later
Y_train <- cbind(yl, yu)
rm(yl, yu)
X_train <- model.matrix(~., data = train_dat[, 3:77])

fit_our <- fsnet(Y = log(Y_train),
                 X = X_train,
                 lam = exp(seq(0, -10, length.out = 11)),
                 alpha = 1,
                 pen_factor = c(rep(0, 7), rep(1, ncol(X_train) - 7)),
                 distr = "ee",
                 fix_var = TRUE,
                 maxit = rep(1e4, 3),
                 tol = rep(1e-7, 2),
                 nfold = 5)

fit_our_small <- fsnet(Y = log(Y_train),
                       X = X_train[, 1:7],
                       lam = 0,
                       distr = "ee",
                       fix_var = TRUE,
                       maxit = rep(1e4, 3),
                       tol = rep(1e-7, 2))

# Calculate mis-classification rate
X_test <- model.matrix(~., data = test_dat[, 3:77])
yl <- 3 * floor(test_dat$time / 3) # Lower endpoint of 3-month interval
yu <- yl + 3 # Upper endpoint of interval
yu[test_dat$event == 0] <- Inf # If censored, may have occurred later
Y_test <- log(cbind(yl, yu))
rm(yl, yu)
lp_our <- X_test %*% fit_our$b_star
lp_our_small <- X_test[, 1:7] %*% fit_our_small$beta
mcr_our <- mean((lp_our < Y_test[, 1]) | (lp_our > Y_test[, 2]))
mcr_our_small <- mean((lp_our_small < Y_test[, 1]) | (lp_our_small > Y_test[, 2]))

pdf("~/GitHub/finite-suppl/breast cancer/fig_trace.pdf", width = 8, height = 4)
par(mfrow = c(1, 1))
par(cex.axis = 1.3, cex.lab = 1.3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))

matplot(x = log(fit_our$full_fit$lam), y = t(fit_our$full_fit$beta[-c(1:7), ]), type = "l",
        xlab = expression(log(lambda)), ylab = expression(alpha), lwd = 1.5)
abline(v = -5)
dev.off()



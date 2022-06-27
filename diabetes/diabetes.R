set.seed(33)
library("glmnetcr")
library(icnet)
data("diabetes")

dim(diabetes)

Y <- matrix(0, nrow = nrow(diabetes), ncol = 2)

# Define cut-off points
table(diabetes[, 1])
t1 <- qnorm(8/24)
t2 <- qnorm(15 / 24)

for(ii in 1:nrow(Y)){
  if(diabetes[ii, 1] == "control"){
    Y[ii, ] <- c(-Inf, t1)
  } else if (diabetes[ii, 1] == "impaired fasting glucose"){
    Y[ii, ] <- c(t1, t2)
  } else{
    Y[ii, ] <- c(t2, Inf)
  }
}



# Marginal model
X <- matrix(1, ncol = 1, nrow = 24)
fit_marg <- icnet(Y = Y,
                   X = X,
                   lam = 0,
                   distr = "norm",
                   alpha = 1,
                   fix_var = TRUE,
                   nfold = 1,
                   tol = c(1e-8, 1e-8),
                  maxit = c(1e3, 1e3, 1e3))
# Probability of "control"
b0 <- (t1 - X[1, ] %*% fit_marg[1, 2:(ncol(X) +1)]) / fit_marg[1, 1]
p0 <- pnorm(b0)

# Probability of impaired fasting glucose
a1 <- (t1 - X[10, ] %*% fit_marg[1, 2:(ncol(X) + 1)]) / fit_marg[1, 1]
b1 <- (t2 - X[10, ] %*% fit_marg[1, 2:(ncol(X) + 1)]) / fit_marg[1, 1]
p1 <- pnorm(b1) - pnorm(a1)

# Probability of typ 2 diabetes
a2 <- (t2 - X[24, ] %*% fit_marg[1, 2:(ncol(X) + 1)]) / fit_marg[1, 1]
p2 <- pnorm(a2, lower = F)

c(p0, p1, p2)



# Effect of predictors
lam_seq <- 2^seq(-7, 0)

X <- scale(as.matrix(diabetes[, -1]))

fit_icnet <- icnet(Y = Y,
                   X = X,
                   lam = lam_seq,
                   distr = "norm",
                   alpha = 1,
                   fix_var = TRUE,
                   nfold = 11,
                   tol = c(1e-8, 1e-8))

matplot(x = rev(lam_seq), fit_icnet$full_fit[, 2:(ncol(X) + 1)], type = "l",
        ylab = "Coefficient", xlab = expression(lambda))
abline(v = fit_icnet$lam_star)
#abline(v = lam_seq[3], lty = 2)
#abline(v = lam_seq[6], lty = 2)

order(fit_icnet$b_star, decreasing = TRUE)[1:20]

sum(fit_icnet$b_star != 0)

# Check with probit
Y_prob <- Y
Y_prob[1:8, 2] <- 0
Y_prob[9:24, 1] <- 0
Y_prob[9:24, 2] <- Inf
fit_probit <- icnet(Y = Y_prob,
                  X = matrix(1, ncol = 1, nrow = nrow(Y)),
                  lam = 0,
                  distr = "norm",
                  alpha = 1,
                  fix_var = TRUE,
                  nfold = 1,
                  tol = c(1e-8, 1e-8),
                  maxit = c(1e3, 1e3, 1e3))
glm(Y_prob[, 2] != 0 ~ 1, family = binomial(link = "probit"))
    
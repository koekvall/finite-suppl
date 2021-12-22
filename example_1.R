library(doParallel)
library(doRNG)
library(icnet)
library(latex2exp)
library(survival)
set.seed(35)
cl <- makeCluster(10)
registerDoParallel(cl)

lambdas <- 10^seq(-1.5, -0.1, length.out = 10)
betas <- log(lambdas)
n <- 200
X <- matrix(1, nrow = n, ncol = 1)
reps <- 1e4
out_mat <- matrix(0, nrow = length(lambdas), ncol = 9)
for(jj in 1:length(lambdas)){
  res_mat <- foreach(ii = 1:reps, .combine = rbind,
                     .errorhandling = "remove", .packages = c("icnet", "survival")) %dorng% {
                       # Censored from above
                       Y <- icnet::generate_ee(X = X, b = betas[jj], d = 1, ymax = 10)
                       beta_hat <- icnet::icnet(Y = Y[, 1:2], X = X, b = betas[jj], lam = 0, alpha = 1, distr = "ee")[1, 1]
                       lambda_hat <- exp(beta_hat)
                       lambda_tilde <- 1 / exp(coef(survreg(Surv(pmin(Y[, 1] + 1, 10), Y[, 1] < 10) ~ 1, dist = "exponential"))[1])
                       res <- c(lambda_hat, lambda_tilde)
                       
                       # Effectively not censored from above
                       Y <- icnet::generate_ee(X = X, b = betas[jj], d = 1, ymax = 1e5)
                       beta_hat <- icnet::icnet(Y = Y[, 1:2], X = X, b = betas[jj], lam = 0, alpha = 1, distr = "ee")[1, 1]
                       lambda_hat <- exp(beta_hat)
                       lambda_tilde <- 1 / exp(coef(survreg(Surv(pmin(Y[, 1] + 1, 10), Y[, 1] < 10) ~ 1, dist = "exponential"))[1])
                       c(res, lambda_hat, lambda_tilde)
                     }
  out_mat[jj, ] <- c(colMeans(res_mat) - lambdas[jj], sqrt(n) * apply(res_mat, 2, sd), lambdas[jj])
  cat("Completed ", jj, "of ", length(lambdas), "\n")
}
colnames(out_mat) <- c("bias", "bias_m", "bias_nc", "bias_m_nc", "sd", "sd_m", "sd_nc", "sd_mnc", "lambda")
stopCluster(cl)

# Plot
PDF <- TRUE
if(PDF) pdf("~/Dropbox/Apps/Overleaf/hd_finite/fig_ex1.pdf", width = 10, height = 4)
par(mfrow = c(1, 2))
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1), oma = c(0, 0, 1, 0))

# Correct
plot(y = out_mat[, 1], x = log10(out_mat[, 9]), type = "l", lwd = 2, ylim = c(-0.3, 1),
     ylab = "", xlab = TeX('log$_{10}(\\lambda)$'), cex.lab = 1.5, cex.axis = 1.5)
points(y = out_mat[, 1], x = log10(out_mat[, 9]), pch = 1)

lines(y = out_mat[, 5], x = log10(out_mat[, 9]), lwd = 2, lty = 1)
points(y = out_mat[, 5], x = log10(out_mat[, 9]), pch = 2)

# Continuous
lines(y = out_mat[, 2], x = log10(out_mat[, 9]), type = "l", lwd = 2, lty = 2)
points(y = out_mat[, 2], x = log10(out_mat[, 9]), pch = 1)
lines(y = out_mat[, 6], x = log10(out_mat[, 9]), lwd = 2, lty = 2)
points(y = out_mat[, 6], x = log10(out_mat[, 9]), pch = 2)

legend("topleft", c("correct", "incorrect", "bias", expression(paste("s.d. times ", sqrt(n)))),
       pch = c(NA, NA, 1, 2), lty = c(1, 2, NA, NA), lwd = c(2, 2, 1, 1), cex = 1.3, bty = "n")

# RMSE
plot(x = log10(out_mat[, 9]), y = sqrt(out_mat[, 2]^2 + out_mat[, 6]^2/n) /
       sqrt(out_mat[, 1]^2 + out_mat[, 5]^2/n),
     type = "l", lwd = 2, ylab = "RMSE ratio", xlab = TeX('log$_{10}(\\lambda)$'), cex.lab = 1.5, cex.axis = 1.5)

# # Correct
# plot(y = out_mat[, 3], x = log10(out_mat[, 9]), type = "l", lwd = 2, ylim = c(-0.3, 1),
#      ylab = "", xlab = TeX('log$_{10}(\\lambda)$'))
# points(y = out_mat[, 3], x = log10(out_mat[, 9]), pch = 1)
# 
# lines(y = out_mat[, 7], x = log10(out_mat[, 9]), lwd = 2, lty = 1)
# points(y = out_mat[, 7], x = log10(out_mat[, 9]), pch = 2)
# 
# # Continuous
# lines(y = out_mat[, 4], x = log10(out_mat[, 9]), type = "l", lwd = 2, lty = 2)
# points(y = out_mat[, 4], x = log10(out_mat[, 9]), pch = 1)
# lines(y = out_mat[, 8], x = log10(out_mat[, 9]), lwd = 2, lty = 2)
# points(y = out_mat[, 8], x = log10(out_mat[, 9]), pch = 2)
# 
# legend("topleft", c("correct", "continuous", "bias", expression(paste("s.d. times ", sqrt(n)))),
#        pch = c(NA, NA, 1, 2), lty = c(1, 2, NA, NA), lwd = c(2, 2, 1, 1))
# 
# plot(x = log10(out_mat[, 9]), y = sqrt(out_mat[, 4]^2 + out_mat[, 8]^2/n) /
#        sqrt(out_mat[, 3]^2 + out_mat[, 7]^2/n),
#      type = "l", lwd = 2, ylab = "RMSE ratio", xlab = TeX('log$_{10}(\\lambda)$'))
# lines(x = log10(out_mat[, 9]), y = sqrt(out_mat[, 4]^2 + out_mat[, 8]^2/n),
#       type = "l", lwd = 2, lty = 2)
# mtext("With censoring", side = 3, line = -1, outer=TRUE, cex=1.2, font = 2)
# mtext("Without censoring", side= 1, line = -17.5, outer=TRUE, cex=1.2, font = 2)
if(PDF) dev.off()

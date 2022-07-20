set.seed(33)
library("glmnetcr")
library(fsnet)
data("diabetes")

dim(diabetes)

# Model effect of predictors
lam_seq <- 2^seq(0, -10)

X <- scale(as.matrix(diabetes[, -1]))
tictoc::tic()
fit_fsnet <- fsnet_cat(Y = diabetes$y,
                   X = X,
                   lam = lam_seq,
                   distr = "norm",
                   alpha = 1,
                   nfold = 10,
                   tol = c(1e-6, 1e-6))
tictoc::toc()

tictoc::tic()
fsnet_cat(Y = diabetes$y,
                       X = X,
                       lam = lam_seq,
                       distr = "norm",
                       alpha = 1,
                       nfold = 1,
                       tol = c(1e-6, 1e-6))
tictoc::toc()


# Non-zero coefficients
num_active <- sum(fit_fsnet$beta_star != 0)
order_active <- order(abs(fit_fsnet$beta_star), decreasing = TRUE)[1:num_active]
fit_fsnet$beta_star[order_active]
colnames(diabetes[, -1])[order_active]
pdf("~/GitHub/finite-suppl/diabetes/fig_trace_diabet.pdf", width = 8, height = 4)
par(mfrow = c(1, 1))
par(cex.axis = 1.3, cex.lab = 1.3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
matplot(x = log(lam_seq), t(fit_fsnet$full_fit$beta), type = "l",
        ylab = "Coefficient", xlab = expression(log(lambda)), lwd = 1.5)
abline(v = log2(fit_fsnet$lam_star))
dev.off()


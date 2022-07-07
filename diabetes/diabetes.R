set.seed(33)
library("glmnetcr")
library(fsnet)
data("diabetes")

dim(diabetes)


# Marginal model is equivalent to unconstrained categorical
fit_marg <- fsnet_cat(Y = diabetes$y, lam = 0, distr = "norm")
# Probability of "control"
p0 <- pnorm(fit_marg$gam[1])

# Probability of impaired fasting glucose
p1 <- pnorm(sum(fit_marg$gam[1:2])) - pnorm(fit_marg$gam[1])

# Probability of typ 2 diabetes
p2 <- pnorm(sum(fit_marg$gam[1:2]), lower = F)

c(p0, p1, p2)
prop.table(table(diabetes$y))

# Model effect of predictors
lam_seq <- 2^seq(0, -10)

X <- scale(as.matrix(diabetes[, -1]))

fit_fsnet <- fsnet_cat(Y = diabetes$y,
                   X = X,
                   lam = lam_seq,
                   distr = "norm",
                   alpha = 1,
                   nfold = 10,
                   tol = c(1e-6, 1e-6))

# Non-zero coefficients
sum(fit_fsnet$beta_star != 0)
order(abs(fit_fsnet$beta_star), decreasing = TRUE)[1:4]

pdf("~/GitHub/finite-suppl/diabetes/fig_trace_diabet.pdf", width = 8, height = 4)
par(mfrow = c(1, 1))
par(cex.axis = 1.3, cex.lab = 1.3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
matplot(x = log(lam_seq), t(fit_fsnet$full_fit$beta), type = "l",
        ylab = "Coefficient", xlab = expression(log(lambda)), lwd = 1.5)
abline(v = log(fit_fsnet$lam_star))
dev.off()


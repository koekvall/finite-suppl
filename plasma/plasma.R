library(fsnet)
plasma_dat <- readRDS("/Users/karekv/GDrive/Research/finite/Data/plasma_dat.Rds")
n <- nrow(plasma_dat)

###############################################################################
# Exponential model
###############################################################################
Y <- matrix(0, nrow = n, ncol = 2)
for(ii in 1:n){
  level_ii = which(plasma_dat$lpa_cat[ii] == levels(plasma_dat$lpa_cat))
  Y[ii, ] <- c((level_ii - 1) * 10, level_ii * 10)
}
Y[Y[, 1] == 120, 2] <- Inf
fit_norm <- fsnet_ic(Y = log(Y), X = matrix(1, nrow = n, ncol = 1), lam = 0,
                   fix_var = FALSE,
                   distr = "norm")

fit_ee <- fsnet_ic(Y = log(Y), X = matrix(1, nrow = n, ncol = 1), lam = 0,
                        fix_var = FALSE,
                        distr = "ee")

bic_norm <- -2 * fit_norm$loglik + log(n) * length(fit_norm$theta)
bic_ee <- -2 * fit_ee$loglik + log(n) * length(fit_ee$theta)

###############################################################################
# Compare to saturated categorical
###############################################################################

fit_cat <- fsnet_cat(Y = plasma_dat$lpa_cat,
                                      distr = "norm",
                                      lam = 0)

bic_cat <- -2 * fit_cat$loglik + log(n) * length(fit_cat$theta)

c("norm" = bic_norm, "ee" = bic_ee, "cat" = bic_cat)

###############################################################################
# Plot estimates of category probabilities
###############################################################################
pdf("~/GitHub/finite-suppl/plasma/fig_pdf_plasma.pdf", width = 8, height = 4)
props <- rbind(prop.table(table(plasma_dat$lpa_cat)),
               plnorm(c(seq(10, 120, by = 10), Inf), meanlog = fit_norm$beta,
                      sdlog = fit_norm$sigma) -
                      plnorm(c(seq(0, 120, by = 10)),
                             meanlog = fit_norm$beta, sdlog = fit_norm$sigma))
rownames(props) = c("Sample proportions", "Censored regression")

barplot(props, ylim = c(0, 0.35),
        ylab = "Estimated probability", xlab = "Lp(a)", beside = TRUE,
        legend = TRUE)
dev.off()
###############################################################################
# Model with predictors
###############################################################################

X <- model.matrix(lpa_cat ~ age * (. - age), data = plasma_dat)

fit_ic <- fsnet_ic(Y = log(Y),
                   X = X,
                   lam = 0,
                   fix_var = FALSE,
                   distr = "norm")

# Test interaction using LRT
fit_noint <- fsnet_ic(Y = log(Y),
                      X = X[, 1:7],
                      lam = 0,
                      fix_var = FALSE,
                      distr = "norm")

pchisq(2 * (fit_ic$loglik - fit_noint$loglik),
       df = length(fit_ic$theta) - length(fit_noint$theta),
       lower = F)

# Coefficients and SEs
est <- c(fit_noint$sigma, fit_noint$beta)
se <- sqrt(diag(solve(-fsnet::hessian_fsnet(Y = log(Y),
                                            X = X[, 1:7],
                                            s = fit_noint$sigma,
                                            b = fit_noint$beta,
                                            fix_s = FALSE,
                                            distr = "norm"))))
pvals <- pchisq((est - c(1, rep(0, length(est) - 1)))^2 / se^2, df = 1, lower = F)

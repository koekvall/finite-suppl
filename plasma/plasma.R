library(icnet)
library(splines2)
plasma_dat <- readRDS("/Users/karekv/GDrive/Research/finite/Data/plasma_dat.Rds")
n <- nrow(plasma_dat)

###############################################################################
# Saturated marginal models
###############################################################################

# Fit both models

fit_marg_ee <- icnet_cat(Y = plasma_dat$lpa_cat,
                                      distr = "ee",
                                      lam = 0)

fit_marg_norm <- icnet_cat(Y = plasma_dat$lpa_cat,
                             distr = "norm",
                             lam = 0)

# Saturated models give the same likelihood regardless of latent density
fit_marg_ee$loglik
fit_marg_norm$loglik

# Saturated models are equivalent to categorical model with maximum likelihood:

# Probability of [0, 10]
pnorm(fit_marg_norm$theta[1])
1 - exp(-exp(fit_marg_ee$theta[1]))
prop.table(table(plasma_dat$lpa_cat))[1] # MLE in categorical model is sample prop.

# Probability of (10, 30]
pnorm(sum(fit_marg_norm$theta[1:2])) - pnorm(sum(fit_marg_norm$theta[1]))
exp(-exp(sum(fit_marg_ee$theta[1]))) - exp(-exp(sum(fit_marg_ee$theta[1:2])))
prop.table(table(plasma_dat$lpa_cat))[2]

# Probability of (120, Inf]
pnorm(sum(fit_marg_norm$theta[1:3]), lower = F)
exp(-exp(sum(fit_marg_ee$theta[1:3])))
prop.table(table(plasma_dat$lpa_cat))[4]

###############################################################################
# Model with predictors
###############################################################################

# Remove intercept, otherwise Z does not have full rank
X <- model.matrix(~ gender * age, data = plasma_dat)[, -1]

fit_ee <- icnet_cat(Y = plasma_dat$lpa_cat,
                        X = X,
                        distr = "ee",
                        lam = 0)

fit_norm <- icnet_cat(Y = plasma_dat$lpa_cat,
                          X = X,
                          distr = "norm",
                          lam = 0)

# Check which model is the better fit
# Note: Models are not nested and have same number of parameters)
fit_ee$loglik
fit_norm$loglik


# Test interaction using LRT
fit_norm_noint <- icnet_cat(Y = plasma_dat$lpa_cat,
                                X = X[, -3],
                                b = fit_norm_new$beta[-3, ],
                                gam = fit_norm_new$gam,
                                distr = "norm",
                                lam = 0)

1 - pchisq(2 * (fit_norm$loglik - fit_norm_noint$loglik), df = 1)

# Test gender using LRT
fit_norm_nogen <- icnet_cat(Y = plasma_dat$lpa_cat,
                                X = X[, 2, drop = F],
                                b = fit_norm_new$beta[-c(1, 3), ],
                                gam = fit_norm_new$gam,
                                distr = "norm",
                                lam = 0)

1 - pchisq(2 * (fit_norm$loglik - fit_norm_nogen$loglik), df = 2)

# Compare to Fisher's exact test (which does not account for age)
table(plasma_dat$lpa_cat, plasma_dat$gender)
fisher.test(table(plasma_dat$lpa_cat, plasma_dat$gender))

# Test age using LRT
fit_norm_noage <- icnet_cat(Y = plasma_dat$lpa_cat,
                                X = X[, 1, drop = F],
                                b = fit_norm_new$beta[-c(1, 2), ],
                                gam = fit_norm_new$gam,
                                distr = "norm",
                                lam = 0)

1 - pchisq(2 * (fit_norm$loglik - fit_norm_noage$loglik), df = 2)

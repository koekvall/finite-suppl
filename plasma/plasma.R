library(icnet)
library(splines2)
plasma_dat <- readRDS("/Users/karekv/GDrive/Research/finite/Data/plasma_dat.Rds")
n <- nrow(plasma_dat)

###############################################################################
# Saturated marginal models
###############################################################################
Z <- matrix(0, nrow = 2 * n, ncol = 3) # Model matrix
M <- matrix(0, nrow = n, ncol = 2) # Offset matrix
for(ii in 1:n){
  if(plasma_dat$lpa_cat[ii] == "[0,10]"){
    M[ii, 1] <- -Inf
    Z[2 * (ii - 1) + 2, ] <- c(1, 0, 0)
  } else if(plasma_dat$lpa_cat[ii] == "(10,30]"){
    Z[2 * (ii - 1) + 1, ] <- c(1, 0, 0)
    Z[2 * (ii - 1) + 2, ] <- c(1, 1, 0)
  } else if(plasma_dat$lpa_cat[ii] == "(30,120]"){
    Z[2 * (ii - 1) + 1, ] <- c(1, 1, 0)
    Z[2 * (ii - 1) + 2, ] <- c(1, 1, 1)
  } else if(plasma_dat$lpa_cat[ii] == "(120,Inf]"){
    M[ii, 2] <- Inf
    Z[2 * (ii - 1) + 1, ] <- c(1, 1, 1)
  } else{
    warning("Missing category")
  }
}

# Constrain positive parameters to ensure a < b
box_constr <- cbind(rep(-Inf, ncol(Z)), rep(Inf, ncol(Z)))
box_constr[2:ncol(Z), 1] <- 0

# Fit both models
fit_marg_ee <- icnet_flex(M = M,
                          Z = Z,
                          theta = c(0, rep(0.1, ncol(Z) - 1)),
                          box_constr = box_constr,
                          verbose = FALSE,
                          method = "prox_newt",
                          distr = "ee",
                          lam = 0)

fit_marg_norm <- icnet_flex(M = M,
                          Z = Z,
                          theta = c(0, rep(0.1, ncol(Z) - 1)),
                          box_constr = box_constr,
                          verbose = FALSE,
                          method = "prox_newt",
                          distr = "norm",
                          lam = 0)

fit_marg_ee_new <- icnet_cat(Y = plasma_dat$lpa_cat,
                                      distr = "ee",
                                      lam = 0)

fit_marg_norm_new <- icnet_cat(Y = plasma_dat$lpa_cat,
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
X <- kronecker(model.matrix(~ gender * age, data = plasma_dat)[, -1],
               -c(1, 1))
# No predictors for infinite ends
X[2 * which(!is.finite(M[, 1])) - 1, ] <- 0
X[2 * which(!is.finite(M[, 2])), ] <- 0

Z_X <- cbind(Z, X)
box_constr_full <- rbind(box_constr, matrix(c(-Inf, Inf), nrow = ncol(X), ncol = 2,
                                            byrow = T))

fit_ee <- icnet_flex(M = M, Z = Z_X,
                       theta = c(fit_marg_ee$theta[1:ncol(Z)], rep(0, ncol(X))),
                       box_constr = box_constr_full,
                       method = "prox_newt",
                       distr = "ee",
                       lam = 0)

fit_ee_new <- icnet_cat(Y = plasma_dat$lpa_cat,
                        X = model.matrix(~ gender * age, data = plasma_dat)[, -1],
                        distr = "ee",
                        lam = 0)

fit_norm <-  icnet_flex(M = M, Z = Z_X,
                        theta = c(fit_marg_norm$theta[1:ncol(Z)], rep(0, ncol(X))),
                        box_constr = box_constr_full,
                        method = "prox_newt",
                        distr = "norm",
                        lam = 0)

fit_norm_new <- icnet_cat(Y = plasma_dat$lpa_cat,
                          X = model.matrix(~ gender * age, data = plasma_dat)[, -1],
                          distr = "norm",
                          lam = 0)

# Check which model is the better fit
# Note: Models are not nested and have same number of parameters)
fit_ee$loglik
fit_norm$loglik


# Test interaction using LRT
fit_norm_noint <- icnet_flex(M = M, Z = Z_X[, -6],
                             theta = fit_norm$theta[1:5],
                             box_constr = box_constr_full[-6, ],
                             method = "prox_newt",
                             distr = "norm",
                             lam = 0)

fit_norm_new_noint <- icnet_cat(Y = plasma_dat$lpa_cat,
                                X = model.matrix(~ gender * age, data = plasma_dat)[, -c(1, 4)],
                                b = fit_norm_new$beta[-3, ],
                                gam = fit_norm_new$gam,
                                distr = "norm",
                                lam = 0)

1 - pchisq(2 * (fit_norm$loglik - fit_norm_noint$loglik), df = 1)
1 - pchisq(2 * (fit_norm_new$loglik - fit_norm_new_noint$loglik), df = 1)

# Test gender using LRT
fit_norm_nogen <- icnet_flex(M = M, Z = Z_X[, -c(4, 6)],
                             theta = fit_norm$theta[c(1:3, 5)],
                             box_constr = box_constr_full[-c(4, 6), ],
                             method = "prox_newt",
                             distr = "norm",
                             lam = 0)
fit_norm_new_nogen <- icnet_cat(Y = plasma_dat$lpa_cat,
                                X = model.matrix(~ gender * age,
                                                 data = plasma_dat)[, -c(1, 2, 4),
                                                                    drop = FALSE],
                                b = fit_norm_new$beta[-c(1, 3), ],
                                gam = fit_norm_new$gam,
                                distr = "norm",
                                lam = 0)

1 - pchisq(2 * (fit_norm$loglik - fit_norm_nogen$loglik), df = 2)
1 - pchisq(2 * (fit_norm_new$loglik - fit_norm_new_nogen$loglik), df = 2)

# Compare to Fisher's exact test (which does not account for age)
table(plasma_dat$lpa_cat, plasma_dat$gender)
fisher.test(table(plasma_dat$lpa_cat, plasma_dat$gender))

# Test age using LRT
fit_norm_noage <- icnet_flex(M = M, Z = Z_X[, -c(5, 6)],
                             theta = fit_norm$theta[1:4],
                             box_constr = box_constr_full[-c(5, 6), ],
                             method = "prox_newt",
                             distr = "norm",
                             lam = 0)
fit_norm_new_noage <- icnet_cat(Y = plasma_dat$lpa_cat,
                                X = model.matrix(~ gender * age,
                                                 data = plasma_dat)[, -c(1, 3, 4),
                                                                    drop = FALSE],
                                b = fit_norm_new$beta[-c(1, 2), ],
                                gam = fit_norm_new$gam,
                                distr = "norm",
                                lam = 0)

1 - pchisq(2 * (fit_norm$loglik - fit_norm_noage$loglik), df = 2)
1 - pchisq(2 * (fit_norm_new$loglik - fit_norm_new_noage$loglik), df = 2)

library(icnet)
library(survival)
n <- 500
X <- matrix(1, nrow = n, ncol = 1)
b_star <- 0.5
Y <- icnet::generate_ee(X = X, b = b_star, d = 1, ymax = 1)

fit_weib <- survreg(Surv(Y[, 1] + 1) ~ 1) 
fit_exp <- survreg(Surv(Y[, 1] + 1) ~ 1, dist = "exponential")

ll_diff <- fit_weib$loglik[1] - fit_exp$loglik[1]

fit_our_exp <- icnet(Y = log(Y[, 1:2]), X = X, lam = 0, fix_s = T,
                 distr = "ee",
                 verbose = F,
                 method = "fista",
                 b = 1,
                 s = 1,
                 L = 10,
                 maxit = 1e4)

# Fit our model with Weibull dist (i.e., coefficient for log(time))
fit_our_weib <- icnet(Y = log(Y[, 1:2]), X = X, lam = 0, fix_s = F,
                 distr = "ee",
                 verbose = F,
                 method = "fista",
                 b = fit_our[1, 2],
                 s = 1,
                 L = 10,
                 maxit = 1e4)


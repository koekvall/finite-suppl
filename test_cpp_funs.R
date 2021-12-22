set.seed(333)
###############################################################################
# Test function of ab
###############################################################################
extreme_cdf <- function(w){1 - exp(-exp(w))}
loglik_ee <- function(ab){log(extreme_cdf(ab[2]) - extreme_cdf(ab[1]))}
loglik_norm <- function(ab){log(pnorm(ab[2]) - pnorm(ab[1]))}
ab_test <- sort(rnorm(2, sd = 1))

# Value
c("naive_ee" = loglik_ee(ab_test),
  "cpp_ee" =  icnet:::loglik_ab(ab_test[1], ab_test[2], 1, 1)[1])

c("naive_norm" = loglik_norm(ab_test),
  "cpp_norm" =  icnet:::loglik_ab(ab_test[1], ab_test[2], 1, 2)[1])

# Gradient
rbind("naive_ee" = numDeriv::grad(loglik_ee, ab_test),
  "cpp_ee" =  icnet:::loglik_ab(ab_test[1], ab_test[2], 2, 1)[2:3])

rbind("naive_norm" = numDeriv::grad(loglik_norm, ab_test),
      "cpp_norm" =  icnet:::loglik_ab(ab_test[1], ab_test[2], 2, 2)[2:3])

# Hessian
rbind("naive_ee" = c(numDeriv::hessian(loglik_ee, ab_test))[-3],
      "cpp_ee" =  icnet:::loglik_ab(ab_test[1], ab_test[2], 2, 1)[4:6])

rbind("naive_norm" = c(numDeriv::hessian(loglik_norm, ab_test))[-3],
      "cpp_norm" =  icnet:::loglik_ab(ab_test[1], ab_test[2], 2, 2)[4:6])
###############################################################################


###############################################################################
# Test objective function
###############################################################################
p <- 40
n <- 20
X <- matrix(rnorm(n * p), ncol = p)
Z <- kronecker(-X, c(1, 1))
b <- matrix(0, nrow = p)
b_test <- runif(p, -0.1, 0.1)
s_test <- runif(1, 0.5, 1)
theta_test <-  b_test / s_test
lam_test <- 0.5
alpha_test <- 0.5
pen_test <- runif(p)

Y_ee <- generate_ee(X = X, b = b_test)
Y_norm <- generate_norm(X = X, b = b_test)

# Value
# ee
ab_ee <- icnet:::get_ab(Z = Z, theta = theta_test, M = log(Y_ee[, 1:2]))
val_cpp_ee <- icnet:::obj_fun(a = ab_ee[, 1], b = ab_ee[, 2],
                              theta = theta_test,
                lam1 = pen_test * alpha_test * lam_test,
                lam2 = pen_test * (1 - alpha_test) * lam_test, dist = 1)

val_cpp_ee_2 <- icnet:::obj_diff_cpp(Z = Z,
                     theta = theta_test,
                     M = log(Y_ee[, 1:2]),
                     lam1 = alpha_test *  lam_test * pen_test,
                     lam2 = (1 - alpha_test) * lam_test * pen_test,
                     order = 0,
                     dist = 1)$obj
naive_obj_ee <- function(theta){
  ab <- icnet:::get_ab(Z = Z, theta = theta, M = log(Y_ee[, 1:2]))
  -mean(log(extreme_cdf(ab[, 2]) - extreme_cdf(ab[, 1]))) +
    lam_test * alpha_test * sum(abs(pen_test * theta)) + 
    lam_test * (1 - alpha_test) * 0.5 * sum(pen_test * theta^2)
}
val_naive_ee <- naive_obj_ee(theta_test)
c("naive_ee" = val_naive_ee, "cpp_ee" = val_cpp_ee, "cpp_ee_2" = val_cpp_ee_2)

# norm
ab_norm <- icnet:::get_ab(Z = Z, theta = theta_test, M = Y_norm[, 1:2])
val_cpp_norm <- icnet:::obj_fun(a = ab_norm[, 1], b = ab_norm[, 2],
                              theta = theta_test,
                              lam1 = pen_test * alpha_test * lam_test,
                              lam2 = pen_test * (1 - alpha_test) * lam_test, dist = 2)

val_cpp_norm_2 <- icnet:::obj_diff_cpp(Z = Z,
                                     theta = theta_test,
                                     M = Y_norm[, 1:2],
                                     lam1 = alpha_test *  lam_test * pen_test,
                                     lam2 = (1 - alpha_test) * lam_test * pen_test,
                                     order = 0,
                                     dist = 2)$obj
naive_obj_norm <- function(theta){
  ab <- icnet:::get_ab(Z = Z, theta = theta, M = Y_norm[, 1:2])
  -mean(log(pnorm(ab[, 2]) - pnorm(ab[, 1]))) +
    lam_test * alpha_test * sum(abs(pen_test * theta)) + 
    lam_test * (1 - alpha_test) * 0.5 * sum(pen_test * theta^2)
}
val_naive_norm <- naive_obj_norm(theta_test)
c("naive_norm" = val_naive_norm, "cpp_norm" = val_cpp_norm, "cpp_norm_2" = val_cpp_norm_2)

# Loglik gradient
naive_ll_ee <- function(theta){
  ab_ee <- icnet:::get_ab(Z = Z, theta = theta, M = log(Y_ee[, 1:2]))
  sum(log(extreme_cdf(ab_ee[, 2]) - extreme_cdf(ab_ee[, 1])))
}
ab_diffs_ee <- icnet:::loglik_ab(ab_ee[, 1], ab_ee[, 2], 1, 1)
rbind("naive_ee" = numDeriv::grad(naive_ll_ee, theta_test),
      "cpp_ee" = c(icnet:::loglik_grad(Z, ab_diffs_ee)))

# Sub-gradient
rbind("naive_ee" = numDeriv::grad(naive_obj_ee, theta_test),
      "cpp_ee" = c(icnet:::obj_diff_cpp(Z = Z,
                                        theta = theta_test,
                                      M = log(Y_ee[, 1:2]),
                                      lam1 = alpha_test *  lam_test * pen_test,
                                      lam2 = (1 - alpha_test) * lam_test * pen_test,
                                      order = 1,
                                      dist = 1)$sub_grad))

rbind("naive_norm" = numDeriv::grad(naive_obj_norm, theta_test),
      "cpp_norm" = c(icnet:::obj_diff_cpp(Z = Z,
                                        theta = theta_test,
                                        M = Y_norm[, 1:2],
                                        lam1 = alpha_test *  lam_test * pen_test,
                                        lam2 = (1 - alpha_test) * lam_test * pen_test,
                                        order = 1,
                                        dist = 2)$sub_grad))

# Hessian
max(abs(icnet:::obj_diff_cpp(Z = Z,
                     theta = theta_test,
                     M = log(Y_ee[, 1:2]),
                     lam1 = alpha_test *  lam_test * pen_test,
                     lam2 = (1 - alpha_test) * lam_test * pen_test,
                     order = 2,
                     dist = 1)$hess -
  numDeriv::hessian(naive_obj_ee, theta_test)))

max(abs(icnet:::obj_diff_cpp(Z = Z,
                             theta = theta_test,
                             M = Y_norm[, 1:2],
                             lam1 = alpha_test *  lam_test * pen_test,
                             lam2 = (1 - alpha_test) * lam_test * pen_test,
                             order = 2,
                             dist = 2)$hess -
          numDeriv::hessian(naive_obj_norm, theta_test)))

###############################################################################
# Test constrained L1 solver
###############################################################################
a <- runif(1); b <- runif(1)
q_fun <- function(x){a * x + 0.5 * b * x^2}
plot(q_fun, -10, 10)

# Unconstrained problem with solution at zero
obj_one <- function(x, lam){q_fun(x) + lam * abs(x)}
plot(function(x)obj_one(x, lam = 0.5), -10, 10)
icnet:::solve_constr_l1(a = a, b = b, c1 = -Inf, c2 = Inf, lam = 0.5)

# Constrained problem with solution at zero
icnet:::solve_constr_l1(a = a, b = b, c1 = 0.1, c2 = Inf, lam = 0.5)
icnet:::solve_constr_l1(a = a, b = b, c1 = -10, c2 = -0.2, lam = 0.5)

# Unconstrained problem with non-zero solution
plot(function(x)obj_one(x, lam = 0.05), -10, 10)
optimize(function(x)obj_one(x, lam = 0.05), interval = c(-10, 10))
icnet:::solve_constr_l1(a = a, b = b, c1 = -Inf, c2 = Inf, lam = 0.05)

# Constrained problem with non-zero solution
plot(function(x)obj_one(x, lam = 0.05), -10, 10)
icnet:::solve_constr_l1(a = a, b = b, c1 = -5, c2 = -4, lam = 0.05)

# Constrained problem with non-zero, positive solution
q_fun <- function(x){-a * x + 0.5 * b * x^2}
obj_one <- function(x, lam){q_fun(x) + lam * abs(x)}
plot(function(x)obj_one(x, lam = 0.05), -10, 10)
icnet:::solve_constr_l1(a = -a, b = b, c1 = 3, c2 = 5, lam = 0.05)


###############################################################################
# Test quadratic approximation
###############################################################################
q_approx_naive_ee <- function(theta, theta_k){
  H <- numDeriv::hessian(naive_ll_ee, theta_k)
  g <- numDeriv::grad(naive_ll_ee, theta_k)
  val <- naive_ll_ee(theta_k)
  theta <- theta - theta_k
  return(val+ crossprod(g, theta) + 0.5 * crossprod(theta, H %*% theta))
}
theta_old <- theta_test + rnorm(p, sd = 0.1)
ab_old_ee <- icnet:::get_ab(Z = Z, theta = theta_old, M = log(Y_ee[, 1:2]))
ab_diffs_old_ee <- icnet:::loglik_ab(ab_old_ee[, 1], ab_old_ee[, 2], 2, 1)
linpred <- icnet:::get_eta(Z, theta_test)
linpred_old <- icnet:::get_eta(Z, theta_old)
icnet:::quad_appr_ll(linpred = linpred, linpred_old = linpred_old, ab_diffs = ab_diffs_old_ee)
q_approx_naive_ee(theta_test, theta_old)


###############################################################################
# Test Newton step
###############################################################################
alpha_newt <- 0
lam_newt <- 0.5 ## Needed when p > n

theta_new_cpp <- icnet:::newton_step(Z = Z,
                    ab = ab_old_ee,
                    ab_diffs = ab_diffs_old_ee,
                    lam1 = rep(0, p),
                    lam2 = rep(lam_newt, p),
                    theta = theta_old,
                    constr = cbind(rep(-Inf, p), rep(Inf, p)),
                    maxit = 100,
                    tol = 1e-8,
                    verbose = F,
                    dist = 1)

theta_new_naive <- theta_old - solve(numDeriv::hessian(naive_ll_ee, theta_old) - diag(n * lam_newt, p),
                                     numDeriv::grad(naive_ll_ee, theta_old) - n * lam_newt * theta_old)

-q_approx_naive_ee(theta_new_cpp, theta_old) + n * 0.5 * lam_newt * sum(theta_new_cpp^2)
-q_approx_naive_ee(theta_new_naive, theta_old) + n * 0.5 * lam_newt * sum(theta_new_naive^2)

###############################################################################
# Test Proximal Newton and FISTA
###############################################################################
library(icnet)
p <- 200
n <- 100
X <- matrix(rnorm(n * p), ncol = p)
Z <- kronecker(-X, c(1, 1))
b <- matrix(0, nrow = p)
b_test <- runif(p, -0.1, 0.1)
s_test <- runif(1, 0.5, 1)
theta_test <-  b_test / s_test
alpha_test <- 0.5
pen_test <- rep(1, p)

Y_ee <- generate_ee(X = X, b = b_test)
Y_norm <- generate_norm(X = X, b = b_test)
lam_seq <- 2^seq(2, -10)

fit_fista_ee <- icnet::icnet(Y = log(Y_ee[, 1:2]),
                             X = X,
                             lam = lam_seq,
                             pen_factor = pen_test,
                             alpha = alpha_test,
                             method = "fista",
                             maxit = 1e4,
                             tol = 1e-6,
                             nfold = 5,
                             distr = "ee")

fit_prox_ee <- icnet::icnet(Y = log(Y_ee[, 1:2]),
                            X = X, lam = lam_seq,
                            pen_factor = pen_test,
                            alpha = alpha_test,
                            method = "prox_newt",
                            nfold = 5,
                            tol = rep(1e-8, 2),
                            distr = "ee",
                            verbose = F,
                            maxit = rep(100, 3))

lam_idx <- 11

obj_prox_ee <- icnet::obj_fun_icnet(Y = log(Y_ee[, 1:2]),
                     X = X,
                     lam = lam_seq[lam_idx],
                     alpha = alpha_test,
                     pen_factor = c(0, pen_test),
                     b = fit_prox_ee$full_fit[lam_idx, 2:(p + 1)],
                     s = 1,
                     distr = "ee",
                     order = 0)$obj

obj_fista_ee <- icnet::obj_fun_icnet(Y = log(Y_ee[, 1:2]),
                     X = X,
                     lam = lam_seq[lam_idx],
                     alpha = alpha_test,
                     pen_factor = c(0, pen_test),
                     b = fit_fista_ee$full_fit[lam_idx, 2:(p + 1)],
                     s = 1,
                     distr = "ee",
                     order = 0)$obj


fit_fista_norm <- icnet::icnet(Y = Y_norm[, 1:2],
                               X = X,
                               lam = lam_seq,
                               pen_factor = pen_test,
                               alpha = alpha_test,
                               method = "fista",
                               maxit = 1e4,
                               tol = 1e-6,
                               nfold = 5,
                               distr = "norm")

fit_prox_norm <- icnet::icnet(Y = Y_norm[, 1:2],
                              X = X,
                              lam = lam_seq,
                              pen_factor = pen_test,
                              alpha = alpha_test,
                              method = "prox_newt",
                              nfold = 5,
                              tol = rep(1e-8, 2))

obj_prox_norm <- icnet::obj_fun_icnet(Y = Y_norm[, 1:2],
                     X = X,
                     lam = lam_seq[lam_idx],
                     alpha = alpha_test,
                     pen_factor = c(0, pen_test),
                     b = fit_prox_norm$full_fit[lam_idx, 2:(p + 1)],
                     s = 1,
                     distr = "norm",
                     order = 0)$obj
  
obj_fista_norm <- icnet::obj_fun_icnet(Y = Y_norm[,1:2],
                     X = X,
                     lam = lam_seq[lam_idx],
                     alpha = alpha_test,
                     pen_factor = c(0, pen_test),
                     b = fit_fista_norm$full_fit[lam_idx, 2:(p + 1)],
                     s = 1,
                     distr = "norm",
                     order = 0)$obj

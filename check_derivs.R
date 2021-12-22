# Calculates difference between numderiv gradient/hessian and explicit gradient/hessian expression of the likelihood function
set.seed(3)
library(numDeriv)
library(icnet)
# Generate data
p <- 4
n <- 20
X <- matrix(rnorm(n * p), ncol = p)
b <- matrix(0, nrow = p)

# Test derivatives for these settings
b_test <- runif(p, -0.1, 0.1)
lam_test <- 0.5
alpha_test <- 0.5
pen_test <- runif(p)

###############################################################################
# Latent normal
###############################################################################
Y <- generate_norm(X, b, d = 0.25, ymax = 0.5, ymin = -0.25)

# Test derivatives w.r.t. end points
loglik <- function(ab){log(pnorm(ab[2]) - pnorm(ab[1]))}
ab_test <- c(-0.1, 0.3)
cbind(c(loglik(ab_test), numDeriv::grad(loglik, ab_test), c(numDeriv::hessian(loglik, ab_test)[c(1, 2, 4)])),
c(icnet:::loglik_ab(ab_test[1], ab_test[2], 2, 2)))

# Manual derivative calculations
naive_obj <- function(b){
  -mean(log(pnorm(Y[, 2] - X %*% b) - pnorm(Y[, 1] - X %*% b))) +
    lam_test * alpha_test * sum(abs(pen_test * b)) + 
    lam_test * (1 - alpha_test) * 0.5 * sum(pen_test * b^2)
}

my_obj <- function(b) icnet:::obj_diff_cpp(Z = kronecker(-X, c(1, 1)),
                               theta = b,
                               M = Y[, 1:2],
                               lam1 = alpha_test *  lam_test * pen_test,
                               lam2 = (1 - alpha_test) * lam_test * pen_test,
                               order = 0, 2)$obj

# Test derivatives with respect to beta
cat("Difference in obj to naive calculation is",  
    my_obj(b_test) - naive_obj(b_test),
    "\n")


cat("Largest difference between numerical and analytical gradient is", 
    max(abs(numDeriv::grad(my_obj, b_test) -
              icnet:::obj_diff_cpp(Z = kronecker(-X, c(1, 1)),
                             theta = b_test,
                             M = Y[, 1:2],
                             lam1 = alpha_test *  lam_test * pen_test,
                             lam2 = (1 - alpha_test) * lam_test * pen_test,
                             order = 1,
                             dist = 2)$sub_grad
          )),
    "\n")

cat("Largest difference between numerical and analytical Hessian is", 
    max(abs(numDeriv::hessian(my_obj, b_test) - 
              icnet:::obj_diff_cpp(Z = kronecker(-X, c(1, 1)),
                                   theta = b_test,
                                   M = Y[, 1:2],
                                   lam1 = alpha_test *  lam_test * pen_test,
                                   lam2 = (1 - alpha_test) * lam_test * pen_test,
                                   order = 2, 2)$hessian)),
    "\n")
###############################################################################


###############################################################################
# Latent Extreme value (~Exponential with log-link)
###############################################################################
Y <- generate_ee(X , b, ymax = 2)

extr_cdf <- function(w){1 - exp(-exp(w))}
# Test derivatives w.r.t. end points
loglik <- function(ab){log(extr_cdf(ab[2]) - extr_cdf(ab[1]))}
ab_test <- c(-0.1, 0.1)

# Objective
loglik(ab_test)
icnet:::loglik_ab(ab_test[1], ab_test[2], 1, 1)[1]

# Gradient w.r.t a, b
numDeriv::grad(loglik, ab_test)
icnet:::loglik_ab(ab_test[1], ab_test[2], 1, 1)[2:3]

# Hessian w.r.t. a, b
numDeriv::hessian(loglik, ab_test)[c(1, 2, 4)]
icnet:::loglik_ab(ab_test[1], ab_test[2], 2, 2)[4:6]


my_obj <- function(theta)icnet:::obj_diff_cpp(Z = kronecker(-X, c(1, 1)),
                                                  theta = theta,
                                                  M = log(Y[, 1:2]),
                                                  lam1 = alpha_test *  lam_test * pen_test,
                                                  lam2 = (1 - alpha_test) * lam_test * pen_test,
                                                  order = 0,
                                           dist = 1)$obj
cat("Difference in obj to naive calculation is",  
    my_obj(b_test) - 
      (-mean(log(pexp(Y[, 2],
                      rate = exp(-X %*% b_test)) - pexp(Y[, 1],rate = exp(-X %*% b_test)))) +
         lam_test * alpha_test * sum(abs(pen_test * b_test)) + 
         lam_test * (1 - alpha_test) * 0.5 * sum(pen_test * b_test^2)),
    "\n")

cat("Largest difference between numerical and analytical gradient is", 
    max(abs(numDeriv::grad(my_obj, b_test) - 
              icnet:::obj_diff_cpp(Z = kronecker(-X, c(1, 1)),
                                   theta = b_test,
                                   M = log(Y[, 1:2]),
                                   lam1 = alpha_test *  lam_test * pen_test,
                                   lam2 = (1 - alpha_test) * lam_test * pen_test,
                                   order = 1,
                                   dist = 1)$sub_grad)),
    "\n")

cat("Largest difference between numerical and analytical Hessian is", 
    max(abs(numDeriv::hessian(my_obj, b_test) - 
              icnet:::obj_diff_cpp(Z = kronecker(-X, c(1, 1)),
                                   theta = b_test,
                                   M = log(Y[, 1:2]),
                                   lam1 = alpha_test *  lam_test * pen_test,
                                   lam2 = (1 - alpha_test) * lam_test * pen_test,
                                   order = 2,
                                   dist = 1)$hessian)),
    "\n")
###############################################################################

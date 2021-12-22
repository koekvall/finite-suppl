library(icnet)
set.seed(44)

# MSE in Exponential example  
n <- 50
b0 <- 1
X <- matrix(1, n, 1)

nreps <- 1000
res_mat <- matrix(0, nreps, 2)
colnames(res_mat) <- c("icnet", "glm")
for(ii in 1:nreps){
  Y <- icnet::generate_ee(X = X, b = b0, d = 1, ymax = 10)
  y_glm <- Y[, 1] + 1
  fit_glm <-  glm(y_glm ~ 1, family = Gamma(link = log))
  fit_icnet <- icnet(Y = log(Y[, 1:2]), X = X, lam = 0, distr = "ee")
  res_mat[ii, ] <- c(fit_icnet[1, 2], coef(fit_glm))
}

# Bias
colMeans(res_mat) - b0

# Standard deviation
sqrt(diag(cov(res_mat)))

# RMSE
sqrt(colMeans((res_mat - b0)^2))


# Mean function in exponential example
beta_seq <- seq(-5, 5, length.out = 10)
true_mean <- beta_seq

for(ii in 1:length(beta_seq)){
  dat <- generate_ee(X = matrix(1, nrow = 1e6, ncol = 1), b = beta_seq[ii], d = 1, ymax = 10)
  
}




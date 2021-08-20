library(icnet)
set.seed(44)
# No predictors
n <- 100
X <- matrix(1, nrow = n, ncol = 1)
b <- 0.1
# Y can take values in +/- (0, 0.5), (0.5, 1), (1, Inf)
# So support has carnality 6
Y <- icnet::generate_norm(X = X, b = b, d = 1, ymax = 0, ymin = 0)
barplot(table(Y[, 1]))

fit <- icnet::icnet(Y = Y[, 1:2], X = X, lam = 0, b = 0.1)

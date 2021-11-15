library(icnet)
set.seed(44)
# No predictors
n <- 100
X <- matrix(1, nrow = n, ncol = 1)
b <- 1
# Y can take values in +/- (0, 0.5), (0.5, 1), (1, Inf)
# So support has carnality 6
Y <- icnet::generate_norm(X = X, b = b, d = 1, ymax = 0, ymin = 0)

barplot(table(Y[, 1]))

fit <- icnet::icnet(Y = Y[, 1:2], X = X, lam = 0, b = 0.1)
glm_fit <- glm(1*(Y[, 2] > 0) ~ 1, family = binomial(link = "probit"))

# Check they are ~ the same
fit[, 1] - coef(glm_fit)


# With predictors
p <- 200
X <- matrix(rnorm(n * p), ncol = p, nrow = n)
X <- scale(X)
b <- c(1, 2, 0.5, rep(0, p - 3))
Y <- icnet::generate_norm(X = X, b = b, d = 1, ymax = 0, ymin = 0)
lam_seq <- 2^seq(5, -5, length.out = 10)
fit <- icnet::icnet(Y = Y[, 1:2], X = X, lam = lam_seq, pen_factor = rep(1, p), alpha = 1)
fit_glmnet <- glmnet(x = X, y = 1*(Y[, 2] > 0), family = binomial(link = "probit"),
                     lambda = lam_seq, intercept = F, alpha = 1, standardize = FALSE)
fit[, 1:p] - t(fit_glmnet$beta)

# Compare objective achieved
print(icnet::obj_diff(y = Y[, 1], X = X, b = fit[10, 1:p], yupp = Y[, 2],
                lam = lam_seq[10], pen_factor = rep(1, p), order = 0,
                dist = "norm")$obj -
icnet::obj_diff(y = Y[, 1], X = X, b = fit_glmnet$beta[, 10], yupp = Y[, 2],
                lam = lam_seq[10], pen_factor = rep(1, p), order = 0,
                dist = "norm")$obj)

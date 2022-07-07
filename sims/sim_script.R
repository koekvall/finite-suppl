named.list <- function(...) {
  l <- list(...)
  names(l) <- as.character( match.call()[-1] )
  l
}

misclass_rate <- function(pred, Y){
  mean((pred < Y[, 1]) | (pred >= Y[, 2]))
}


do_one_sim <- function(set)
{
  set.seed(set$seed)
  R_X <- chol(0.5^abs(outer(1:set$p, 1:set$p, "-")))
  X <- cbind(1, scale(matrix(rnorm(set$n * set$p), nrow = set$n, ncol = set$p) %*% R_X))

  # Three predictors for low-dim exponential
  b0 <- c(1, 0.5, -0.5)
  Y <- fsnet::generate_ee(X[, 1:3], b = b0, d = set$d, ymax = set$ymax)
  y_glm <- Y[, 1] + set$d
  fit_exp_our <- fsnet::fsnet(Y = log(Y[, 1:2]),
                              X = X[, 1:3],
                              lam = 0,
                              fix_var = T,
                              s = set$s,
                              method = "prox_newt",
                              maxit = rep(100, 3),
                              distr = "ee")
  fit_exp_glm <- glm(y_glm ~ 0 + X[, 1:3], family = Gamma(link = log))
  b_ld <- fit_exp_our$beta
  b_glm <- coef(fit_exp_glm)
  pred_ld <- exp(X[, 1:3] %*% b_ld)
  pred_glm <- exp(X[, 1:3] %*% b_glm)
  mcr_ld <- misclass_rate(pred_ld, Y[, 1:2])
  mcr_glm <- misclass_rate(pred_glm, Y[, 1:2])



  # p predictors but no intercept for high-dim
  b0_hd <- c(b0, rep(0, set$p - 3))
  Y <- fsnet::generate_norm(X = X[, -1, drop = F],
                            b = b0_hd,
                            d = set$d,
                            ymax = set$ymax,
                            ymin = set$ymin,
                            sigma = set$s)
  y_glm <- Y[, 1]

  fit_glmnet <- glmnet::cv.glmnet(x = X[, -1, drop = F],
                         y = y_glm,
                         lambda = set$lam_seq,
                         alpha = 1,
                         intercept = FALSE,
                         standardize = FALSE,
                         nfolds = 5)
  lam_glmnet <- fit_glmnet$lambda.min

  fit_norm_prox <- fsnet::fsnet(Y = Y[, 1:2],
                               X = X[, -1, drop = F],
                               lam = set$lam_seq,
                               alpha = 1,
                               b = b0_hd,
                               s = set$s,
                               fix_var = TRUE,
                               method = "prox_newt",
                               distr = "norm",
                               tol = c(1e-8, 1e-8),
                               verbose = F,
                               nfold = 5)
  lam_our <- fit_norm_prox$lam_star



  b_hd <- fit_norm_prox$b_star
  pred_hd <- X[, -1] %*% b_hd
  pred_glmnet <- predict(fit_glmnet, s = "lambda.min", newx = X[, -1])
  b_glmnet <- coef(fit_glmnet, s = "lambda.min")[-1]
  Y_new <- fsnet::generate_norm(X = X[, -1, drop = F], b = b0_hd, d = set$d,
                            ymax = 100, ymin = -100, sigma = 1)
  mcr_hd <- misclass_rate(pred_hd, Y_new[, 1:2])
  mcr_glmnet <- misclass_rate(pred_glmnet, Y_new[, 1:2])
  return(named.list(b_ld, b_glm, b0, b_hd, b_glmnet, b0_hd, mcr_ld, mcr_glm, mcr_hd, mcr_glmnet, lam_our,
                    lam_glmnet))
}

###############################################################################
# Simulation
###############################################################################
library(doParallel)
library(doRNG)
today <- as.numeric(format(Sys.time(), "%H%d%m%y"))
seed_start <- 10 * today
cl <- makeCluster(10)
registerDoParallel(cl)
n_sims <- 500

base_set <- list(n = 100, p = 200, d = 1, ymin = -5, ymax = 5,
                 lam_seq = 2^seq(2, -15, length.out = 10), seed = 35, s = 1)
d_list <- seq(0.5, 5, by = 0.5)

out_list <- list()
idx <- 1
for(ii in 1:length(d_list)){
  out_list[[ii]] <- foreach(kk = 1:n_sims,
                            .combine = rbind,
                            .packages = c("glmnet", "fsnet"),
                            .errorhandling = "remove") %dorng%{
                              set_kk <- base_set
                              set_kk$seed <- seed_start + kk
                              set_kk$d <- d_list[ii]
                              c(do_one_sim(set_kk), set_kk)
                            }
  cat("Completed simulation ", idx, "of ", length(d_list), "\n")
  idx <- idx + 1
  seed_start <- seed_start + n_sims
}
stopCluster(cl)
res_mat <- do.call(rbind, out_list)
saveRDS(res_mat, paste0("~/GitHub/finite-suppl/sims/", today, ".Rds"))

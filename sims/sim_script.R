named.list <- function(...) {
  l <- list(...)
  names(l) <- as.character( match.call()[-1] )
  l
}

# set <- list(n = 50, p = 100, d = 1, ymin = -5, ymax = 5, lam_seq = 2^seq(2, -15, length.out = 10), seed = 35, s = 0.5)

do_one_sim <- function(set)
{
  set.seed(set$seed)
  R_X <- chol(0.5^abs(outer(1:set$p, 1:set$p, "-")))
  X <- cbind(1, scale(matrix(rnorm(set$n * set$p), nrow = set$n, ncol = set$p) %*% R_X))

  # Three predictors for low-dim exponential
  b0 <- c(1, 0.5, -0.5)
  Y <- icnet::generate_ee(X[, 1:3], b = b0, d = set$d, ymax = set$ymax)
  y_glm <- Y[, 1] + set$d
  fit_exp_our <- icnet::icnet(Y = log(Y[, 1:2]),
                              X = X[, 1:3],
                              lam = 0,
                              fix_var = F,
                              method = "prox_newt",
                              maxit = rep(100, 3),
                              distr = "ee")
  fit_exp_glm <- glm(y_glm ~ 0 + X[, 1:3], family = Gamma(link = log))
  b_our_ld <- fit_exp_our[1, 2:4]
  b_glm <- coef(fit_exp_glm)

  # p predictors but no intercept for high-dim
  b0_hd <- c(b0, rep(0, set$p - 3))
  Y <- icnet::generate_norm(X = X[, -1, drop = F],
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

  fit_norm_prox <- icnet::icnet(Y = Y[, 1:2],
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


  misclass_rate <- function(eta, Y){
      mean((eta < Y[, 1]) | (eta >= Y[, 2]))
  }

  b_our_hd <- fit_norm_prox$b_star
  pred <- X[, -1] %*% b_our_hd
  pred_glmnet <- predict(fit_glmnet, s = "lambda.min", newx = X[, -1])
  b_glmnet <- coef(fit_glmnet, s = "lambda.min")[-1]
  Y_new <- icnet::generate_norm(X = X[, -1, drop = F], b = b0_hd, d = set$d,
                            ymax = 100, ymin = -100, sigma = 1)
  mcr <- misclass_rate(pred, Y_new[, 1:2])
  mcr_glmnet <- misclass_rate(pred_glmnet, Y_new[, 1:2])
  return(named.list(b_our_ld, b_glm, b0, b_our_hd, b_glmnet, b0_hd, mcr, mcr_glmnet, lam_our = fit_norm_prox$lam_star,
                    lam_glmnet = fit_glmnet$lambda.min))
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
n_sims <- 100

base_set <- list(n = 50, p = 100, d = 1, ymin = -5, ymax = 5,
                 lam_seq = 2^seq(2, -15, length.out = 10), seed = 35, s = 1)
d_list <- seq(0.5, 5, by = 0.5)

out_list <- list()
idx <- 1
for(ii in 1:length(d_list)){
  out_list[[ii]] <- foreach(kk = 1:n_sims,
                            .combine = rbind,
                            .packages = c("glmnet", "icnet"),
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

###############################################################################
# Figures
###############################################################################
out_mat <- readRDS("~/GitHub/finite-suppl/sims/15211221.Rds")
res_mat <- matrix(0, nrow = nrow(out_list), 19)
colnames(res_mat) <- c( "ssb_our", "sse_our", "ssb_glm", "sse_glm", "ssb_hd",
"sse_hd", "ssb_glmnet", "sse_glmnet", "mcr_hd", "mcr_glmnet", "n", "p" , "d",
"ymin", "ymax", "s", "lam_our", "lam_glmnet", "seed")
for(ii in 1:nrow(out_list)){
  res <- out_mat[ii, ]
  res_mat[ii, ] <- c("ssb_our" = sum(res$b_our_ld^2),
                     "sse_our" = sum((res$b_our_ld - res$b0)^2),
                     "ssb_glm" = sum(res$b_glm^2),
                     "sse_glm" = sum((res$b_glm - res$b0)^2),
                     "ssb_hd" = sum(res$b_our_hd^2),
                     "sse_hd" = sum((res$b_our_hd - res$b0_hd)^2),
                     "ssb_glmnet" = sum(res$b_glmnet^2),
                     "sse_glmnet" = sum((res$b_glmnet - res$b0_hd)^2),
                     "mcr_hd" = res$mcr,
                     "mcr_glmnet" = res$mcr_glmnet,
                     "n" = res$n,
                     "p" = res$p,
                     "d" = res$d,
                     "ymin" = res$ymin,
                     "ymax" = res$ymax,
                     "s" = res$s,
                     "lam_our" = res$lam_our,
                     "lam_glmnet" = res$lam_glmnet,
                     "seed" = res$seed
                     )
}


as_tibble(res_mat) %>% select("sse_our", "sse_glm", "d") %>%
  pivot_longer(cols = c(sse_our, sse_glm),
               values_to = "sse",
               names_prefix = "sse_",
               names_to = "method") %>%
  group_by(d, method) %>%
  summarize(mse = mean(sse),
            lower = mean(sse) - 1.96 * sd(sse) / n(),
            upper = mean(sse) + 1.96 * sd(sse) / n(),
            n = n()) %>%
  ggplot(aes(x = d, y = mse, group = method, col = method)) + geom_point() + geom_line() +
  geom_smooth(aes(ymin = lower, ymax = upper, fill = method, colour = method), stat = "identity")

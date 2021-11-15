# set <- c(n, p, distr, d, ymin, ymax, lam_min, lam_max, nsim, seed, ncore)
do_one_sim <- function(set)
{
  set.seed(set$seed)
  library(doParallel)
  library(doRNG)
  
  penalize <- ifelse(set$lam_min > 0, TRUE, FALSE)
  X <- scale(matrix(runif(set$n * set$p, -1, 1), nrow = set$n, ncol = set$p))
  b0 <-c(0.5, -1, 2)
  if(set$p > 3){
    b0 <- c(b0, rep(0, set$p - 3))
  }else{
    b0 <- b0[1:set$p]
  }
  
  cl <- makeCluster(set$ncore)
  registerDoParallel(cl)
  
  res_mat <- foreach(ii = 1:set$nsim, .combine = rbind,
                     .errorhandling = "remove", .packages = c("glmnet")) %dorng%
    {
      # Function to calculate prediction misclassification rate
      misclass_rate <- function(eta, Y, distr){
        if(distr == "norm"){
          out <- mean(ifelse((eta < Y[, 1]) | (eta >= Y[, 2]), 1, 0))
        } else{
          out <- mean(ifelse((exp(eta) < Y[, 1]) | (exp(eta) >= Y[, 2]), 1, 0))
        }
        out
      }
      
      if(set$distr == "norm"){
        Y_mat <- icnet::generate_norm(X, b0, d = set$d, ymax = set$ymax, ymin = set$ymin)
        Y_new_mat <- icnet::generate_norm(X, b0, d = set$d, ymax = set$ymax, ymin = set$ymin)
      } else{
        Y_mat <- icnet::generate_ee(X, b0, d = set$d, ymax = set$ymax)
        Y_new_mat <- icnet::generate_ee(X, b0, d = set$d, ymax = set$ymax)
      }
      
      if(penalize){
        lam_seq <- 2^seq(log2(set$lam_max), log2(set$lam_min), length.out = 5)
        fit <- icnet::icnet(Y = Y_mat[, 1:2], X = X, b = b0, lam = lam_seq,
                            alpha = 1, pen_factor = rep(1, set$p),
                            distr = "norm", nfold = 10)
        pred<- X %*% fit$b_star
        lam_star <- fit$lam_star
        
        fit_other <- glmnet::cv.glmnet(x = X,
                                       y = ifelse(Y_mat[, 1] >= 0, Y_mat[, 1], Y_mat[, 2]),
                                       lambda = lam_seq,
                                       alpha = 1,
                                       intercept = FALSE,
                                       standardize = FALSE,
                                       family = "gaussian")
        pred_other <- predict(fit_other, s = "lambda.min", newx = X)
        lam_star_other <- fit_other$lambda.min
        
        mcr <- misclass_rate(pred, Y_new_mat[, 1:2], "norm")
        mcr_other <- misclass_rate(pred_other, Y_new_mat[, 1:2], "norm")
        
      } else if (distr == "norm"){
        fit <- icnet::icnet(Y = Y_mat[, 1:2], X = X, b = b0, lam = 0,
                            distr = "norm")
        pred <- X %*% fit[1, 1:p]
        fit_other <- lm(ifelse(Y_mat[, 1] >= 0, Y_mat[, 1], Y_mat[, 2]) ~ 0 + X)
        pred_other <- X %*% coef(fit_other)
        mcr <- misclass_rate(pred, Y_new_mat[, 1:2], "norm")
        mcr_other <- misclass_rate(pred_other, Y_new_mat[, 1:2], "norm")
        lam_star <- NA
        lam_star_other <- NA
      } else{
        fit <- icnet::icnet(Y = Y_mat[, 1:2], X = X, b = b0, lam = 0,
                            distr = "ee")
        pred <- X %*% fit[1, 1:p]
        fit_other <- glm(Y_mat[, 1] ~ 0 + X, family = Gamma(link = "log"))
        pred_other <- -(X %*% coef(fit_other)) # Sign difference
        mcr <- misclass_rate(pred, Y_new_mat[, 1:2], "ee")
        mcr_other <- misclass_rate(pred_other, Y_new_mat[, 1:2], "ee")
        lam_star <- NA
        lam_star_other <- NA
      }
      
      out <- c("mcr" = mcr, "mcr_other" = mcr_other, "lam_star" = lam_star,
               "lam_star_other" = lam_star_other, unlist(set[1:2]),
               unlist(set[4:11]))
  }
  stopCluster(cl)
  return(res_mat)
}

try_set <- list("n" = 100, "p" = 500, "distr" = "norm", "d" = 1, "ymin" = -100,
                "ymax" = 100, "lam_min" = 2^(-10), "lam_max" = 2^5, "nsim" = 50,
                "seed" = 33, "ncore" = 11)
try_res <- do_one_sim(try_set)
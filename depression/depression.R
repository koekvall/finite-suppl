library(tidyverse)
library(haven)
library(icnet)

dep_dat <- left_join(read_xpt("~/GitHub/finite-suppl/depression/DPQ_J.xpt"),
                     read_xpt("~/GitHub/finite-suppl/depression/DEMO_J.xpt"),
                              .id = "SEQN") %>%
  select(SEQN:DPQ100, RIAGENDR, RIDEXAGM) %>%
  na.omit()

# Sum up points and give zero points for refused answers and "don't know"
dep <- rowSums(dep_dat[, 2:10] * (dep_dat[, 2:10] <= 3))

use_dat <- dep_dat %>% transmute(score = dep,
                                 age = RIDEXAGM / 12,
                                 male = RIAGENDR == 1)

#fit_lm <- lm(score ~ ., data = use_dat)

Y <- cbind(use_dat$score, use_dat$score + 1)
Y[Y[, 1] == 27, 2] <- Inf # None with max score

fit_our <- icnet::icnet(Y = log(Y),
                        X = model.matrix(score ~ ., data = use_dat),
                        lam = 0,
                        fix_var = F,
                        maxit = rep(1e3, 3),
                        tol = rep(1e-10, 2),
                        distr = "ee")

fit_our_exp <- icnet::icnet(Y = log(Y),
                            X = model.matrix(score ~ ., data = use_dat),
                            lam = 0,
                            fix_var = T,
                            maxit = rep(1e3, 3),
                            tol = rep(1e-10, 2),
                            distr = "ee")

# Estimates and standard errors
sb_hat <- fit_our[1, 1:4]

se <- sqrt(diag(solve(-icnet::hessian_icnet(Y = log(Y),
                                     X = model.matrix(score ~ ., data = use_dat),
                                     s = sb_hat[1],
                                     b = sb_hat[2:4],
                                     fix_s = FALSE,
                                     distr = "ee"))))
# p-values
unname(2 * pnorm(abs((sb_hat - c(1, 0, 0, 0)) / se), lower = F))

# Plot estimated latent distribution
x_male <- c(1, median(use_dat$age), 1)
x_fem <- c(1, median(use_dat$age), 0)
pdf <- matrix(0, nrow = 28, ncol = 2)
colnames(pdf) <- c("male", "female")
rownames(pdf) <- 0:27
for(ii in 0:27){
  if(ii == 0){
    a_fem <- a_male <- -Inf
  } else{
    a_male <- (log(ii) - crossprod(sb_hat[2:4], x_male)) / sb_hat[1]
    a_fem <- (log(ii) - crossprod(sb_hat[2:4], x_fem)) / sb_hat[1]
  }
  
  if(ii == 27){
    b_male <- b_fem <- Inf
  } else{
    b_male <- (log(ii + 1) - crossprod(sb_hat[2:4], x_male)) / sb_hat[1]
    b_fem <- (log(ii + 1) - crossprod(sb_hat[2:4], x_fem)) / sb_hat[1]
  }
  
  pdf[ii + 1, ] <- exp(c(icnet:::loglik_ab(a = a_male, b = b_male, order = 0, dist = 1),
                     icnet:::loglik_ab(a = a_fem, b = b_fem, order = 0, dist = 1)))
}


# Create plots
# pdf("~/GitHub/finite-suppl/depression/fig_pdf.pdf", width = 8, height = 4)
par(cex.axis = 1.3, cex.lab = 1.3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
barplot(t(pdf), beside = T, legend.text = TRUE, ylab = "Probability", xlab = "Score")
# dev.off()

#pdf("~/GitHub/finite-suppl/depression/fig_hist.pdf", width = 8, height = 4)
par(cex.axis = 1.3, cex.lab = 1.3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
hist(use_dat$score, xlab = "Score", ylab = "Count", main = "")
#dev.off()


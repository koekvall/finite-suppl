###############################################################################
# Figures
###############################################################################
library(tidyverse)
library(ggthemes)
library(cowplot)

out_mat <- readRDS("~/GitHub/finite-suppl/sims/13080722.Rds")
res_mat <- matrix(0, nrow = nrow(out_mat), 17)
colnames(res_mat) <- c("sse_ld", "mcr_ld", "sse_glm",  "mcr_glm", "sse_hd",
                         "mcr_hd", "sse_glmnet", "mcr_glmnet", "n", "p" , "d",
                        "ymin", "ymax", "s", "lam_our", "lam_glmnet", "seed")
for(ii in 1:nrow(out_mat)){
  res <- out_mat[ii, ]
  res_mat[ii, ] <- c("sse_ld" = sum((res$b_ld - res$b0)^2),
                     "mcr_ld" = res$mcr_ld,
                     "sse_glm" = sum((res$b_glm - res$b0)^2),
                     "mcr_glm" = res$mcr_glm,
                     "sse_hd" = sum((res$b_hd[1:3] - res$b0_hd[1:3])^2),
                     "mcr_hd" = res$mcr_hd,
                     "sse_glmnet" = sum((res$b_glmnet[1:3] - res$b0_hd[1:3])^2),
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

# LD MSE
p_ld_mse <- as_tibble(res_mat) %>% select("sse_ld", "sse_glm", "d") %>%
  filter(d  < 5) %>%
  pivot_longer(cols = c(sse_ld, sse_glm),
               values_to = "sse",
               names_prefix = "sse_",
               names_to = "method") %>%
  mutate(method = recode(method,"ld" = "fsnet","glmnet" = "glmnet")) %>%
  group_by(d, method) %>%
  summarize(mse = mean(sse),
            sd_sse = sd(sse),
            n = n()) %>%
  ggplot(aes(x = d, y = mse, group = method, linetype = method)) +
  geom_ribbon(aes(ymin = mse - 1.96 * sd_sse / sqrt(n),
                  ymax = mse + 1.96 * sd_sse / sqrt(n)),
                  fill = "grey70") +
  geom_point() + geom_line() +
  theme_Publication() +
  labs(x = "Interval size", y = "MSE") +
  ggtitle("Low dimension") +
  theme(legend.title = element_blank())

# LD MCR

p_ld_mcr <- as_tibble(res_mat) %>% select("mcr_ld", "mcr_glm", "d") %>%
  filter(d  < 5) %>%
  pivot_longer(cols = c(mcr_ld, mcr_glm),
               values_to = "mcr",
               names_prefix = "mcr_",
               names_to = "method") %>%
  mutate(method = recode(method,"ld" = "fsnet","glmnet" = "glmnet")) %>%
  group_by(d, method) %>%
  summarize(mmcr = mean(mcr),
            sd_mcr = sd(mcr),
            n = n()) %>%
  ggplot(aes(x = d, y = mmcr, group = method, linetype = method)) +
  geom_ribbon(aes(ymin = mmcr - 1.96 * sd_mcr / sqrt(n),
                  ymax =  mmcr + 1.96 * sd_mcr / sqrt(n)),
              fill = "grey70") +
  geom_point() + geom_line() +
  theme_Publication() +
  labs(x = "Interval size", y = "MMCR") +
  ggtitle("Low dimension") +
  theme(legend.title = element_blank())


#HD SSE

p_hd_mse <- as_tibble(res_mat) %>% select("sse_hd", "sse_glmnet", "d") %>%
  filter(d  < 5) %>%
  pivot_longer(cols = c(sse_hd, sse_glmnet),
               values_to = "sse",
               names_prefix = "sse_",
               names_to = "method") %>%
  mutate(method = recode(method,"hd" = "fsnet","glmnet" = "glmnet")) %>%
  group_by(d, method) %>%
  summarize(mse = mean(sse),
            sd_sse = sd(sse), n = n()) %>%
  ggplot(aes(x = d, y = mse, group = method, linetype = method)) +
  geom_ribbon(aes(ymin = mse - 1.96 * sd_sse / sqrt(n),
                  ymax = mse + 1.96 * sd_sse / sqrt(n)),
              fill = "grey70") +
  geom_point() + geom_line() +
  theme_Publication() +
  labs(x = "Interval size", y = "MSE") +
  ggtitle("High dimension") +
  theme(legend.title = element_blank())

#HD MCR

p_hd_mcr <- as_tibble(res_mat) %>% select("mcr_hd", "mcr_glmnet", "d") %>%
  filter(d  < 5) %>%
  pivot_longer(cols = c(mcr_hd, mcr_glmnet),
               values_to = "mcr",
               names_prefix = "mcr_",
               names_to = "method") %>%
  mutate(method = recode(method,"hd" = "fsnet","glmnet" = "glmnet")) %>%
  group_by(d, method) %>%
  summarize(mmcr = mean(mcr),
            sd_mcr = sd(mcr), n = n()) %>%
  ggplot(aes(x = d, y = mmcr, group = method, linetype = method)) +
  geom_ribbon(aes(ymin = mmcr - 1.96 * sd_mcr / sqrt(n),
                  ymax = mmcr + 1.96 * sd_mcr / sqrt(n)),
              fill = "grey70") +
  geom_point() + geom_line() +
  theme_Publication() +
  labs(x = "Interval size", y = "MMCR") +
  ggtitle("High dimension") +
  theme(legend.title = element_blank())




ggsave(filename = "~/Dropbox/Apps/Overleaf/finite_regr/fig_simsR1.pdf",
       plot = plot_grid(p_ld_mse, p_ld_mcr, p_hd_mse, p_hd_mcr, nrow = 2),
       device = "pdf",
       width = 10,
       height = 6,
       units = "in")





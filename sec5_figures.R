library(tidyverse)
library(patchwork)
library(ggthemes)

PDF <- TRUE

# Output and wokring directory
out_dir <- ""
setwd("")
source("ggplot_theme_Publication.R")


###############################################################################
# Main simulation results
###############################################################################

# UNCOMMENT THIS TO READ RAW DATA PRODICED BY SIM
# fig_dat <- list.files(path = "./Data/Raw/Coverage/", pattern = ".Rds", full.names = T) %>%
#   map_dfr(readRDS) %>% pivot_longer(cols = c("LL", "RLL", "LRT", "WLD"),
#                                     names_to = c("stat")) %>%
#   distinct(n1, n2, p, param, seed, stat, type, .keep_all = TRUE) %>%
#   mutate(df = ifelse(type == "corr", 4, 3)) %>%
#   group_by(param, stat, type, n1, n2, p) %>%
#   summarize(prop = mean(value <= qchisq(0.95, df)), reps = n()) %>%
#   mutate(se = sqrt(prop * (1 - prop) / reps))
# saveRDS(object = fig_dat, file = "./Data/fig_dat.Rds")
fig_dat <- readRDS(file = "./Data/fig_dat.Rds")
###############################################################################
# Large n, small p
###############################################################################

p1 <- fig_dat %>% filter(type == "indep", n1 == 1000, p == 2) %>%
  filter(stat %in% c("LRT", "WLD"), param <= 0.2) %>%
  
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9))  + ylim(0.93, 0.985)

p2 <- fig_dat %>% filter(type == "corr", n1 == 1000, p == 2) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Correlation", y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Correlated") + ylim(0.93, 0.97)

p3 <- fig_dat %>% filter(type == "cross",  n1 == 40, p == 2) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed") + ylim(0.87, 1)

p_large_n <- p1 + p2 + p3 + plot_layout(ncol = 3)
if(PDF) pdf(paste0(out_dir, "fig_large_n.pdf"), width = 15, height = 6)
p_large_n
if(PDF) dev.off()


###############################################################################
# Small n, small p
###############################################################################

p4 <- fig_dat %>% filter(type == "indep", param <= 0.5, n1 == 20, n2 == 3, p == 2) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9)) + ggtitle("Independent") +
  geom_hline(yintercept=0.95)

p5 <-  fig_dat %>% filter(type == "corr", param <= -0.5, n1 == 20, n2 == 3, p == 2) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Correlation", y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Correlated")+
  geom_hline(yintercept=0.95)

p6 <- fig_dat %>% filter(type == "cross", param <= 0.5, n1 == 10, p == 2) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed")+
  geom_hline(yintercept=0.95)

p_small_n <- p4 + p5 + p6 + plot_layout(ncol = 3)
if(PDF) pdf(paste0(out_dir, "fig_small_n.pdf"), width = 15, height = 6)
p_small_n
if(PDF) dev.off()

###############################################################################
# Large n, large p
###############################################################################

p7 <- fig_dat %>% filter(type == "indep", param <= 0.15, n1 == 200, p == 400) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9)) + ggtitle("Independent") + ylim(0.94, 0.98)

p8 <- fig_dat %>% filter(type == "corr", param <= -0.85, n1 == 200, p == 400) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Correlation", y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Correlated") + ylim(0.94, 0.98)

p9 <- fig_dat %>% filter(type == "cross", param <= 0.15, n1 == 40, p == 400) %>%
  ggplot(aes(x = param, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed") + ylim(0.94, 0.98)

p_large_n_large_p<- p7 + p8 + p9 + plot_layout(ncol = 3)
if(PDF) pdf(paste0(out_dir, "fig_large_n_large_p.pdf"), width = 15, height = 6)
p_large_n_large_p
if(PDF) dev.off()


# ###############################################################################
# # Small n, large p
# ###############################################################################
#
# p10 <- fig_dat %>% filter(type == "indep", n1 == 20, n2 == 5, p == 20) %>%
#   ggplot(aes(x = param, y = prop, group = stat)) +
#   geom_line(aes(color = stat, lty = stat)) +
#   geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
#   labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "Coverage") +
#   theme_Publication() + theme(legend.position = c(0.9, 0.9)) + ggtitle("Independent") + ylim(0.9, 1)
#
# p11 <- fig_dat %>% filter(type == "corr", n1 == 20, n2 == 5, p == 20) %>%
#   ggplot(aes(x = param, y = prop, group = stat)) +
#   geom_line(aes(color = stat, lty = stat)) +
#   geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
#   labs(lty = "Statistic", color = "Statistic", x = "Correlation", y = "") +
#   theme_Publication() + theme(legend.position = "none") + ggtitle("Correlated") + ylim(0.9, 1)
#
# p12 <- fig_dat %>% filter(type == "cross", param <= 0.15, n1 == 10, p == 20) %>%
#   ggplot(aes(x = param, y = prop, group = stat)) +
#   geom_line(aes(color = stat, lty = stat)) +
#   geom_ribbon(aes(ymin = prop - 2 * se, ymax = prop + 2 * se), alpha = 0.1) +
#   labs(lty = "Statistic", color = "Statistic", x = "Variance", y = "") +
#   theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed") + ylim(0.94, 0.98)
#
# p_small_n_large_p<- p10 + p11 + p12 + plot_layout(ncol = 3)
# if(PDF) pdf(paste0(out_dir, "fig_small_n_large_p.pdf"), width = 15, height = 6)
# p_small_n_large_p
# if(PDF) dev.off()

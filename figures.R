library(tidyverse)
library(patchwork)
library(ggthemes)

PDF <- TRUE

# Output and wokring directory
out_dir <- "~/GitHub/uniform_lmm_suppl/Figures/"
source("~/GitHub/uniform_lmm_suppl/R1/ggplot_theme_Publication.R")


###############################################################################
# Main simulation results
###############################################################################

# UNCOMMENT THIS TO READ RAW DATA PRODICED BY SIMS
fig_dat <- list.files(path = "~/GitHub/uniform_lmm_suppl/Results/",
                      pattern = ".Rds", full.names = T) %>%
  map_dfr(readRDS) %>% pivot_longer(cols = c("RSCR", "PSCR", "LRT", "WLD"),
                                    names_to = c("stat")) %>%
  #distinct(n1, n2, p, param, seed, stat, type, .keep_all = TRUE) %>%
  mutate(df = ifelse(type == "corr", 4, 3)) %>%
  group_by(type, n1, n2, p, psi1, psi2, psi3, psi4, stat) %>%
  summarize(prop = mean(value <= qchisq(0.95, df)), reps = n()) %>%
  mutate(se = sqrt(prop * (1 - prop) / reps))
saveRDS(object = fig_dat, file = "~/GitHub/uniform_lmm_suppl/fig_dat.Rds")
#fig_dat <- readRDS(file = "./Data/fig_dat.Rds")

###############################################################################
# Clustered data, random effect correlation near negative unity
###############################################################################
p1 <- fig_dat %>% filter(type == "corr", n1 == 500, p == 2, psi1 == psi3) %>%
  ggplot(aes(x = psi2, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[2]), y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9)) + ggtitle("Large sample") + ylim(0.925, 0.975)+
  geom_hline(yintercept=0.95)

p2 <-  fig_dat %>% filter(type == "corr", n1 == 20, n2 == 3, p == 2) %>%
  ggplot(aes(x = psi2, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[2]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Small sample")+
  geom_hline(yintercept=0.95)

p3 <- fig_dat %>% filter(type == "corr", n1 == 200, p == 100) %>%
  ggplot(aes(x = psi2, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[2]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Many predictors") + ylim(0.925, 0.975) +
  geom_hline(yintercept=0.95)

p_corr <- p1 + p2 + p3 + plot_layout(ncol = 3)

if(PDF) ggsave(filename = paste0(out_dir, "fig_corr.pdf"), plot = p_corr,
               width = 15, height = 6)

###############################################################################
# Clustered data, random effect variances near zero 
###############################################################################

p4 <- fig_dat %>% filter(type == "indep", n1 == 500, p == 2, psi1 == psi2) %>%
  # filter(stat %in% c("LRT", "WLD"), param <= 0.2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9)) + ggtitle("Large sample")  + 
  ylim(0.93, 0.985) + geom_hline(yintercept=0.95)


p5 <- fig_dat %>% filter(type == "indep", n1 == 20, n2 == 3, p == 2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Small sample") +
  geom_hline(yintercept=0.95)

p6 <- fig_dat %>% filter(type == "indep", n1 == 200, p == 100) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Many predictors") + ylim(0.935, 0.985) +
  geom_hline(yintercept=0.95)

p_indep<- p4 + p5 + p6 + plot_layout(ncol = 3)

if(PDF) ggsave(filename = paste0(out_dir, "fig_indep.pdf"), plot = p_indep,
               width = 15, height = 6)

###############################################################################
# Crossed random effects, random effect variances near zero
###############################################################################

p7 <- fig_dat %>% filter(type == "cross",  n1 == 40, p == 2, psi1 == psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Large sample") + ylim(0.89, 1)+
  geom_hline(yintercept=0.95)


p8 <- fig_dat %>% filter(type == "cross", n1 == 10, p == 2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Small sample") +
  geom_hline(yintercept=0.95)

p9 <- fig_dat %>% filter(type == "cross", n1 == 20, p == 80) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Many predictors") + ylim(0.86, 1) +
  geom_hline(yintercept=0.95)

p_crossed <- p7 + p8 + p9 + plot_layout(ncol = 3)

if(PDF) ggsave(filename = paste0(out_dir, "fig_cross.pdf"), plot = p_crossed,
               width = 15, height = 6)

###############################################################################
# Supplementary simulation results
###############################################################################

###############################################################################
# Different values of psi1 and psi2, independent clusters
###############################################################################

p41 <- p4 + ggtitle("") + ylim(0.93, 0.985)

p42 <- fig_dat %>% filter(type == "indep", n1 == 500, p == 2, psi1 > psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = "none")  + ylim(0.93, 0.985) +
  geom_hline(yintercept=0.95)

p43 <- fig_dat %>% filter(type == "indep", n1 == 500, p == 2, psi1 < psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = "none")  + ylim(0.93, 0.985) +
  geom_hline(yintercept=0.95)

p_diff_psi <- p41 + p42 + p43 + plot_layout(ncol = 3)

if(PDF) ggsave(filename = paste0(out_dir, "fig_diff_psi.pdf"), plot = p_diff_psi,
               width = 15, height = 6)

###############################################################################
# Different values of n1 and n2 with crossed random effects
###############################################################################

p71 <- p7 + ggtitle("") + ylim(0.7, 1) + theme(legend.position = c(0.1, 0.3))

p72 <- fig_dat %>% filter(type == "cross",  n2 == 80, p == 2, psi1 == psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("") + ylim(0.7, 1)+
  geom_hline(yintercept=0.95)

p_diff_ni <- p71 + p72 + plot_layout(ncol = 2)

if(PDF) ggsave(filename = paste0(out_dir, "fig_diff_ni.pdf"), plot = p_diff_ni,
               width = 15, height = 6)

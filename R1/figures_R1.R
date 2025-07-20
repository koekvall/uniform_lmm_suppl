library(tidyverse)
library(patchwork)
library(ggthemes)

PDF <- TRUE

# Output and wokring directory
out_dir <- "~/GitHub/uniform_lmm_suppl/R1/Figures/"
source("~/GitHub/uniform_lmm_suppl/R1/ggplot_theme_Publication.R")


###############################################################################
# Main simulation results
###############################################################################

# UNCOMMENT THIS TO READ RAW DATA PRODICED BY SIMS
fig_dat <- list.files(path = "~/GitHub/uniform_lmm_suppl/R1/Results/",
                      pattern = ".Rds", full.names = T) %>%
  map_dfr(readRDS) %>% pivot_longer(cols = c("RSCR", "PSCR", "LRT", "WLD"),
                                    names_to = c("stat")) %>%
  #distinct(n1, n2, p, param, seed, stat, type, .keep_all = TRUE) %>%
  mutate(df = ifelse(type == "corr", 4, 3)) %>%
  group_by(type, n1, n2, p, psi1, psi2, psi3, psi4, stat) %>%
  summarize(prop = mean(value <= qchisq(0.95, df)), reps = n()) %>%
  mutate(se = sqrt(prop * (1 - prop) / reps))
saveRDS(object = fig_dat, file = "~/GitHub/uniform_lmm_suppl/R1/fig_dat.Rds")
#fig_dat <- readRDS(file = "./Data/fig_dat.Rds")
###############################################################################
# Large n, small p
###############################################################################

p1 <- fig_dat %>% filter(type == "indep", n1 == 1000, p == 2, psi1 == psi2) %>%
  # filter(stat %in% c("LRT", "WLD"), param <= 0.2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9))  + ylim(0.93, 0.985) +
  geom_hline(yintercept=0.95)

p2 <- fig_dat %>% filter(type == "corr", n1 == 1000, p == 2, psi1 == psi3) %>%
  ggplot(aes(x = psi2, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[2]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Correlated") + ylim(0.93, 0.99)+
  geom_hline(yintercept=0.95)

p3 <- fig_dat %>% filter(type == "cross",  n1 == 40, p == 2, psi1 == psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed") + ylim(0.87, 1)+
  geom_hline(yintercept=0.95)

p_large_n <- p1 + p2 + p3 + plot_layout(ncol = 3)
if(PDF) pdf(paste0(out_dir, "fig_large_n.pdf"), width = 15, height = 6)
p_large_n
if(PDF) dev.off()


###############################################################################
# Small n, small p
###############################################################################

p4 <- fig_dat %>% filter(type == "indep", n1 == 20, n2 == 3, p == 2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9)) + ggtitle("Independent") +
  geom_hline(yintercept=0.95)

p5 <-  fig_dat %>% filter(type == "corr", n1 == 20, n2 == 3, p == 2) %>%
  ggplot(aes(x = psi2, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[2]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Correlated")+
  geom_hline(yintercept=0.95)

p6 <- fig_dat %>% filter(type == "cross", n1 == 10, p == 2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed") +
  geom_hline(yintercept=0.95)

p_small_n <- p4 + p5 + p6 + plot_layout(ncol = 3)
if(PDF) pdf(paste0(out_dir, "fig_small_n.pdf"), width = 15, height = 6)
p_small_n
if(PDF) dev.off()

###############################################################################
# Large n, large p
###############################################################################

p7 <- fig_dat %>% filter(type == "indep", n1 == 200, p == 100, stat != "PSCR") %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9)) + ggtitle("Independent") + ylim(0.92, 1) +
  geom_hline(yintercept=0.95)

p8 <- fig_dat %>% filter(type == "corr", n1 == 200, p == 100, stat != "PSCR") %>%
  ggplot(aes(x = psi2, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[2]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Correlated") + ylim(0.90, 0.99) +
  geom_hline(yintercept=0.95)

p9 <- fig_dat %>% filter(type == "cross", n1 == 20, p == 80) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed") + ylim(0.80, 1) +
  geom_hline(yintercept=0.95)

p_large_n_large_p<- p7 + p8 + p9 + plot_layout(ncol = 3)
if(PDF) pdf(paste0(out_dir, "fig_large_n_large_p.pdf"), width = 15, height = 6)
p_large_n_large_p
if(PDF) dev.off()

###############################################################################
# Supplementary simulation results
###############################################################################

###############################################################################
# Different values of psi1 and psi2
###############################################################################
p11 <- fig_dat %>% filter(type == "indep", n1 == 1000, p == 2, psi1 > psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9))  + ylim(0.93, 0.985)+
  geom_hline(yintercept=0.95)

p12 <- fig_dat %>% filter(type == "indep", n1 == 1000, p == 2, psi1 < psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "Coverage") +
  theme_Publication() + theme(legend.position = c(0.9, 0.9))  + ylim(0.93, 0.985)+
  geom_hline(yintercept=0.95)

p_diff_psi <- p1 + p11 + p12 + plot_layout(ncol = 3)
if(PDF) pdf(paste0(out_dir, "fig_diff_psi.pdf"), width = 15, height = 6)
p_diff_psi
if(PDF) dev.off()

###############################################################################
# Different values of n1 and n2 with crossed random effects
###############################################################################

p31 <- fig_dat %>% filter(type == "cross",  n2 == 80, p == 2, psi1 == psi2) %>%
  ggplot(aes(x = psi1, y = prop, group = stat)) +
  geom_line(aes(color = stat, lty = stat)) +
  geom_ribbon(aes(ymin = prop - 2 * se, ymax = pmin(prop + 2 * se, 1)), alpha = 0.1) +
  labs(lty = "Statistic", color = "Statistic", x = expression(psi[1]), y = "") +
  theme_Publication() + theme(legend.position = "none") + ggtitle("Crossed") + ylim(0.7, 1)+
  geom_hline(yintercept=0.95)

p_diff_ni <- p3 + p31 + plot_layout(ncol = 2)
if(PDF) pdf(paste0(out_dir, "fig_diff_ni.pdf"), width = 15, height = 6)
p_diff_ni
if(PDF) dev.off()

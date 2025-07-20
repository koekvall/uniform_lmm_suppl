set.seed(1336)
library(limestest)
library(trust)
library(foreach)
library(doParallel)
library(doRNG)

PDF <- TRUE
out_dir <- "~/GitHub/uniform_lmm_suppl/R1/Figures/"

one_sim <- function(psi) {
  n1 <- 50
  n2 <- 5
  n <- n1 * n2
  # Matrix of covariance params
  Psi1 <- matrix(c(psi[1], psi[2], psi[2], psi[3]), 2, 2)
  Psi <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1)

  ed <- eigen(Psi1)
  R1 <- ed$vectors %*% (sqrt(ed$values) * t(ed$vectors))
  R <- Matrix::kronecker(Matrix::Diagonal(n1), R1)

  H1 <- Matrix::kronecker(Matrix::Diagonal(n1),
                          matrix(c(1, 0, 0, 0), 2, 2))
  H2 <- Matrix::kronecker(Matrix::Diagonal(n1),
                           matrix(c(0, 1, 1, 0), 2, 2))
  H3 <- Matrix::kronecker(Matrix::Diagonal(n1),
                           matrix(c(0, 0, 0, 1), 2, 2))
  H <- cbind(H1, H2, H3)

  # number of covariates (for constructing the matrix X)
  X <- cbind(1, matrix(rbinom(n, 1, 0.5), nrow = n, ncol = 1))


  # Corresponds to random intercept and random slope
  Z <- Matrix::bdiag(lapply(1:n1,
                             function(ii)X[((ii - 1) * n2 + 1):(ii * n2), ,
                                           drop = FALSE]))

  # Random effects model; mean zero ==> e = y
  y <- rnorm(n, sd = sqrt(psi[4])) + Z %*% crossprod(R, rnorm(ncol(R)))

  # Statistic for psi[-1]
  stuff <- limestest:::loglik_psi(Z = Z,
                                 ZtZXe = crossprod(Z, cbind(Z, X, y)),
                                 e = y, H = H, Psi_r = Psi * (1 / psi[4]),
                                 psi_r = psi[4], get_val = TRUE, get_score = TRUE,
                                 get_inf = TRUE, expected = TRUE)
  statboth <- as.numeric(crossprod(stuff$score, solve(stuff$inf_mat,
                                                stuff$score)))
  # Wald statistic
  mc_dat <- data.frame("obs" = as.numeric(y), "v1" = as.numeric(X[, 2]),
                       "clust" = as.factor(rep(1:n1, each = n2)))
  fit <- lme4::lmer(obs ~ 0 + (1 + v1|clust), data = mc_dat)
  psi_hat <- limestest:::get_psi_hat_lmer(fit)
  Psi_hat <- limestest:::Psi_from_H_cpp(psi_hat[-4], H)
  stuff_wald <- limestest:::loglik_psi(Z = Z,
                        ZtZXe = crossprod(Z, cbind(Z, X, y)),
                        e = y, H = H, Psi_r = Psi_hat * (1 / psi_hat[4]),
                        psi_r = psi_hat[4], get_val = TRUE, get_score = TRUE,
                        get_inf = TRUE, expected = FALSE)
  stat_wald <- as.numeric(crossprod(psi_hat - psi,
                                    stuff_wald$inf_mat %*% (psi_hat - psi)))

  stat_lrt <- 2 * (stuff_wald$value - stuff$value)

  c(statboth, stat_wald, stat_lrt)
}

# Settings

# Change fourth element to change error variance
psi_true <- c(1e-3, 0, 1e-3, 1)

nsim <- 10000

### For multi-core simulation
totalCores = detectCores()
# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 1)
registerDoParallel(cluster)

vals <- foreach(r = 1:nsim, .packages = c("trust", "limestest", "Matrix", "lme4"),
                .combine = rbind) %dorng% {
  one_sim(psi_true)
}
stopCluster(cluster)

############ PLOTTING ###############
if(PDF) pdf(paste0(out_dir, "fig_intro.pdf"), width = 15, height = 6)
cex_val <- 1.8
par(mfrow = c(1, 3))
par(mar = c(5, 4.5, 4, 2))
probs <- ppoints(nsim/100)
plot(x = qchisq(probs, df = 4), y = quantile(vals[, 1], probs = probs),
     xlim = c(0, 14), ylim = c(0, 14),
     main = "Score", xlab = "Chi-squared quantile",
     ylab = "Empirical quantile",
     cex.main = cex_val, cex.axis = cex_val,cex.lab = cex_val)
abline(a = 0, b = 1)
plot(x = qchisq(probs, df = 4), y = quantile(vals[, 2], probs = probs),
     xlim = c(0, 14), ylim = c(0, 14),
     main = "Wald", xlab = "Chi-squared quantile",
     ylab = "Empirical quantile",
     cex.main = cex_val, cex.axis = cex_val,cex.lab = cex_val)
abline(a = 0, b = 1)
plot(x = qchisq(probs, df = 4), y = quantile(vals[, 3], probs = probs),
     xlim = c(0, 14), ylim = c(0, 14),
     main = "LRT", xlab = "Chi-squared quantile",
     ylab = "Empirical quantile",
     cex.main = cex_val, cex.axis = cex_val, cex.lab = cex_val)
abline(a = 0, b = 1)
if(PDF) dev.off()
#####################################



cat("Coverage of Score, Wald, and LRT: ", colMeans(vals < qchisq(0.95, df = 4)),
    "\n")




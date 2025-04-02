do_one_indep_sim <- function(inner_seed, n1, n2, X, Z, Psi, psi0, H, R, XtX,
  XtZ, ZtZ, Xb, REML)
{
  # Generate data
  set.seed(inner_seed)
  n <- n1 * n2
  y <- Xb + rnorm(n, sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))

  #############################################################################
  # Score statistics
  #############################################################################
  Psi0 <- Psi / psi0
  XtY <- crossprod(X, y)
  YtZ <- crossprod(y, Z)
  stuff_REML <- limestest::res_ll(XtX =XtX,
                                  XtY = XtY,
                                  XtZ = XtZ,
                                  ZtZ = ZtZ,
                                  YtZ = YtZ,
                                  Y = y,
                                  X = X,
                                  Z = Z,
                                  H = H,
                                  Psi0 = Psi0,
                                  psi0 = psi0,
                                  score = TRUE,
                                  finf = TRUE,
                                  lik = TRUE)
  e <- y - X %*% stuff_REML$beta # REML beta = ML beta at true psi^0
  Zte <- crossprod(Z, e)
  stuff <- limestest::loglik_psi(Z = Z,
                                 ZtZXe = cbind(ZtZ, t(XtZ), Zte),
                                 e = e,
                                 H = H,
                                 Psi0 = Psi0,
                                 psi0 = psi0,
                                 loglik = TRUE,
                                 score = TRUE,
                                 finf = TRUE)
  test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
  test_stat_REML <- as.vector(crossprod(stuff_REML$score,
                                        solve(stuff_REML$finf, stuff_REML$score)))


  #############################################################################
  # Wald statistics. Note: Needs to be before LR statistics to get max loglik
  #############################################################################
  mc_data_frame <- tidyr::tibble("out" = as.vector(y),
                                 "clust" = as.factor(rep(1:n1, each = n2)))
  mc_data_frame <- dplyr::bind_cols(mc_data_frame,
                                    tidyr::as_tibble(X[, -1, drop = F]))
  fit <- lme4::lmer(out ~ . -clust + (1|clust) + (0 + V1|clust),
                    data = mc_data_frame, REML = REML)
  VC <- as.data.frame(lme4::VarCorr(fit))

  # MLEs or REMLEs
  psi_est <- VC$vcov[1:2]
  psi0_est <- VC$vcov[3]
  Psi0_est <- Matrix::Diagonal(2 * n1, x = c(rep(psi_est[1] / psi0_est, n1),
                                             rep(psi_est[2] / psi0_est, n1)))

  # Get observed information (for Wald) and loglik (for LRT) at estimates
  e <- y - X %*% fixef(fit) # Replace residuals
  Zte <- crossprod(Z, e)
  stuff_at_est <- limestest::loglik_psi(Z = Z,
                                ZtZXe = cbind(ZtZ, t(XtZ), Zte),
                                e = e,
                                H = H,
                                Psi0 = Psi0_est,
                                psi0 = psi0_est,
                                loglik = TRUE,
                                score = FALSE,
                                finf = TRUE,
                                expected = FALSE)
  psi <- c(Psi[1, 1], Psi[ncol(Psi), ncol(Psi)])
  est_err<- c(psi0_est, psi_est) - c(psi0, psi)
  test_stat_wald <- as.vector(crossprod(est_err, stuff_at_est$finf %*% est_err))

  #############################################################################
  # LR statistics
  #############################################################################
  if(REML){
    ll_null <- stuff_REML$ll
    ll_max  <- limestest::res_ll(XtX =XtX,
                                  XtY = XtY,
                                  XtZ = XtZ,
                                  ZtZ = ZtZ,
                                  YtZ = YtZ,
                                  Y = y,
                                  X = X,
                                  Z = Z,
                                  H = H,
                                  Psi0 = Psi0_est,
                                  psi0 = psi0_est,
                                  score = FALSE,
                                  finf = FALSE,
                                  lik = TRUE)$ll
  } else{
    ll_null <- stuff$ll
    ll_max  <- stuff_at_est$ll
  }

  test_stat_lrt <- 2 * (ll_max - ll_null)

  # Return
  c(test_stat, test_stat_REML, test_stat_lrt, test_stat_wald, inner_seed)
}




do_one_clust_sim <- function(inner_seed, n1, n2, X, Z, Psi1, psi0, H, R, XtX,
  XtZ, ZtZ, Xb, REML)
{
  # Generate data
  set.seed(inner_seed)
  n <- n1 * n2
  y <- Xb + rnorm(n, sd = sqrt(psi0)) + Z %*% crossprod(R, rnorm(ncol(R)))

  #############################################################################
  # Score statistics
  #############################################################################
  Psi0 <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1 / psi0)
  XtY <- crossprod(X, y)
  YtZ <- crossprod(y, Z)
  stuff_REML <- limestest::res_ll(XtX =XtX,
                                  XtY = XtY,
                                  XtZ = XtZ,
                                  ZtZ = ZtZ,
                                  YtZ = YtZ,
                                  Y = y,
                                  X = X,
                                  Z = Z,
                                  H = H,
                                  Psi0 = Psi0,
                                  psi0 = psi0,
                                  score = TRUE,
                                  finf = TRUE,
                                  lik = TRUE)
  e <- y - X %*% stuff_REML$beta # REML beta = ML beta at true psi^0
  Zte <- crossprod(Z, e)
  stuff <- limestest::loglik_psi(Z = Z,
                                 ZtZXe = cbind(ZtZ, t(XtZ), Zte),
                                 e = e,
                                 H = H,
                                 Psi0 = Psi0,
                                 psi0 = psi0,
                                 loglik = TRUE,
                                 score = TRUE,
                                 finf = TRUE)
  test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
  test_stat_REML <- as.vector(crossprod(stuff_REML$score,
                                        solve(stuff_REML$finf, stuff_REML$score)))
  #############################################################################
  # Wald statistics. Note: Needs to be before LR statistics to get max loglik
  #############################################################################
  mc_data_frame <- tidyr::tibble("out" = as.vector(y),
                                 "clust" = as.factor(rep(1:n1, each = n2)))
  mc_data_frame <- dplyr::bind_cols(mc_data_frame,
                                    tidyr::as_tibble(X[, -1, drop = F]))
  fit <- lme4::lmer(out ~ . - clust + (1 + V1|clust), data = mc_data_frame,
                    REML = REML)

  # MLEs or REMLEs
  VC <- as.data.frame(lme4::VarCorr(fit))
  Psi1_est <- matrix(c(VC$vcov[1], VC$vcov[3], VC$vcov[3], VC$vcov[2]), 2, 2)
  psi0_est <- VC$vcov[4]
  Psi0_est <- Matrix::kronecker(Matrix::Diagonal(n1), Psi1_est / psi0_est)
  psi_est <- c(Psi1_est[1, 1], Psi1_est[2, 1], Psi1_est[2, 2])

  # Get observed information (for Wald) and loglik (for LRT) at estimates
  e <-  y - X %*% fixef(fit) # Replace residuals
  Zte <- crossprod(Z, e)
  stuff_at_est <- limestest::loglik_psi(Z = Z,
                                        ZtZXe = cbind(ZtZ, t(XtZ), Zte),
                                        e = e,
                                        H = H,
                                        Psi0 = Psi0_est,
                                        psi0 = psi0_est,
                                        loglik = TRUE,
                                        score = FALSE,
                                        finf = TRUE,
                                        expected = FALSE)
  psi <- c(Psi1[1, 1], Psi1[2, 1], Psi1[2, 2])
  est_err<- c(psi0_est, psi_est) - c(psi0, psi)
  test_stat_wald <- as.vector(crossprod(est_err, stuff_at_est$finf %*% est_err))

  #############################################################################
  # LR statistics
  #############################################################################
  if(REML){
    ll_null <- stuff_REML$ll
    ll_max  <- limestest::res_ll(XtX =XtX,
                                 XtY = XtY,
                                 XtZ = XtZ,
                                 ZtZ = ZtZ,
                                 YtZ = YtZ,
                                 Y = y,
                                 X = X,
                                 Z = Z,
                                 H = H,
                                 Psi0 = Psi0_est,
                                 psi0 = psi0_est,
                                 score = FALSE,
                                 finf = FALSE,
                                 lik = TRUE)$ll
  } else{
    ll_null <- stuff$ll
    ll_max  <- stuff_at_est$ll
  }

  test_stat_lrt <- 2 * (ll_max - ll_null)

  # Return
  c(test_stat, test_stat_REML, test_stat_lrt, test_stat_wald, inner_seed)
}



do_one_cross_sim <- function(inner_seed, n1, n2, X, Z, psi, psi0, H, XtX, XtZ, ZtZ, Xb, REML)
{
  # Generate data
  set.seed(inner_seed)
  n <- n1 * n2
  Psi <- Matrix::Diagonal(n1 + n2, c(rep(psi[1], n1), rep(psi[2], n2)))
  y <- Xb + rnorm(n, sd = sqrt(psi0)) + Z %*% crossprod(sqrt(Psi), rnorm(ncol(Psi)))

  #############################################################################
  # Score statistics
  #############################################################################
  Psi0 <- Psi / psi0
  XtY <- crossprod(X, y)
  YtZ <- crossprod(y, Z)
  stuff_REML <- limestest::res_ll(XtX =XtX,
                                  XtY = XtY,
                                  XtZ = XtZ,
                                  ZtZ = ZtZ,
                                  YtZ = YtZ,
                                  Y = y,
                                  X = X,
                                  Z = Z,
                                  H = H,
                                  Psi0 = Psi0,
                                  psi0 = psi0,
                                  score = TRUE,
                                  finf = TRUE,
                                  lik = TRUE)
  e <- y - X %*% stuff_REML$beta # REML beta = ML beta at true psi^0
  Zte <- crossprod(Z, e)
  stuff <- limestest::loglik_psi(Z = Z,
                                 ZtZXe = cbind(ZtZ, t(XtZ), Zte),
                                 e = e,
                                 H = H,
                                 Psi0 = Psi0,
                                 psi0 = psi0,
                                 loglik = TRUE,
                                 score = TRUE,
                                 finf = TRUE)
  test_stat <- as.vector(crossprod(stuff$score, solve(stuff$finf, stuff$score)))
  test_stat_REML <- as.vector(crossprod(stuff_REML$score,
                                        solve(stuff_REML$finf, stuff_REML$score)))


  #############################################################################
  # Wald statistics. Note: Needs to be before LR statistics to get max loglik
  #############################################################################
  mc_data_frame <- tidyr::tibble("out" = as.vector(y),
                                 "clust1" = as.factor(rep(1:n1,
                                                          each = n2)),
                                 "clust2" = as.factor(rep(1:n2, n1)))
  mc_data_frame <- dplyr::bind_cols(mc_data_frame,
                                    tidyr::as_tibble(X[, -1, drop = F]))
  fit <- lme4::lmer(out ~ . - clust1 - clust2 + (1|clust1) + (1|clust2),
                    data = mc_data_frame, REML = REML)

  # MLEs or REMLEs
  VC <- as.data.frame(lme4::VarCorr(fit))
  idx <- c(which(VC$grp == "clust1"), which(VC$grp == "clust2"),
           which(VC$grp == "Residual"))
  psi_est <- VC$vcov[idx[1:2]]
  psi0_est <- VC$vcov[idx[3]]
  Psi0_est <- Matrix::Diagonal(n1 + n2, c(rep(psi_est[1] / psi0_est, n1),
                                          rep(psi_est[2] / psi0_est, n2)))

  # Get observed information (for Wald) and loglik (for LRT) at estimates
  e <- y - X %*% fixef(fit) # Replace residuals
  Zte <- crossprod(Z, e)
  stuff_at_est <- limestest::loglik_psi(Z = Z,
                                        ZtZXe = cbind(ZtZ, t(XtZ), Zte),
                                        e = e,
                                        H = H,
                                        Psi0 = Psi0_est,
                                        psi0 = psi0_est,
                                        loglik = TRUE,
                                        score = FALSE,
                                        finf = TRUE,
                                        expected = FALSE)
  est_err<- c(psi0_est, psi_est) - c(psi0, psi)
  test_stat_wald <- as.vector(crossprod(est_err, stuff_at_est$finf %*% est_err))

  #############################################################################
  # LR statistics
  #############################################################################
  if(REML){
    ll_null <- stuff_REML$ll
    ll_max  <- limestest::res_ll(XtX =XtX,
                                 XtY = XtY,
                                 XtZ = XtZ,
                                 ZtZ = ZtZ,
                                 YtZ = YtZ,
                                 Y = y,
                                 X = X,
                                 Z = Z,
                                 H = H,
                                 Psi0 = Psi0_est,
                                 psi0 = psi0_est,
                                 score = FALSE,
                                 finf = FALSE,
                                 lik = TRUE)$ll
  } else{
    ll_null <- stuff$ll
    ll_max  <- stuff_at_est$ll
  }

  test_stat_lrt <- 2 * (ll_max - ll_null)

  # Return
  c(test_stat, test_stat_REML, test_stat_lrt, test_stat_wald, inner_seed)
}

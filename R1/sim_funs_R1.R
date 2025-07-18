do_one_sim <- function(settings)
{
  # Unpack settings
  n1 <- settings$n1
  n2 <- settings$n2
  X <- settings$X
  Z <- settings$Z
  psi <- settings$psi
  H <- settings$H
  Xb <- settings$Xb
  R <- settings$R
  XtX <- settings$XtX
  XtZ = settings$XtZ
  ZtZ <- settings$ZtZ
  REML <- settings$REML
  type = settings$type

  r <- length(psi)
  n <- nrow(Z)

  Psi <- limestest:::Psi_from_H_cpp(psi_mr = psi[-r], H = H)
  q <- ncol(Psi)

  if(is.null(R)){
    R <- chol(Psi)
  }

  if(is.null(XtX)){
    XtX <- crossprod(X)
  }

  if(is.null(XtZ)){
    XtZ <- crossprod(XtZ)
  }

  if(is.null(ZtZ)){
    ZtZ <- crossprod(Z)
    ZtZ <- methods::as(ZtZ, "generalMatrix")
  }
  # Generate data
  y <- as.vector(Xb + rnorm(n, sd = sqrt(psi[r])) + Z %*% crossprod(R, rnorm(q)))

  #############################################################################
  # Score statistics
  #############################################################################
  Psi_r <- Psi / psi[r]
  XtY <- as.matrix(crossprod(X, y))
  ZtY <- as.matrix(crossprod(Z, y))
  
  stuff_REML <- limestest:::res_ll_cpp(Y = y,
                                       X = X,
                                       Z = Z,
                                       XtY = XtY,
                                       ZtY = ZtY,
                                       XtX = XtX,
                                       XtZ = XtZ,
                                       ZtZ = ZtZ,
                                       Psi_r = Psi_r,
                                       psi_r = psi[r],
                                       H = H,
                                       get_val = TRUE,
                                       get_score = TRUE,
                                       get_inf = TRUE)
  e <- as.vector(y - X %*% stuff_REML$beta) # REML beta = ML beta at true psi
  Zte <- as.matrix(crossprod(Z, e))
  stuff <- limestest:::loglik_psi_cpp(e = e, Z = Z, Zte = Zte, XtZ = XtZ,
                                      ZtZ = ZtZ, Psi_r = Psi_r, psi_r = psi[r],
                                      H = H, get_val = TRUE, get_score = TRUE,
                                      get_inf = TRUE, expected = TRUE)
  test_stat <- as.vector(crossprod(stuff$score, solve(stuff$inf_mat, stuff$score)))
  test_stat_REML <- as.vector(crossprod(stuff_REML$score,
                                        solve(stuff_REML$inf_mat, stuff_REML$score)))


  #############################################################################
  # Wald statistics. Note: Needs to be before LR statistics to get max loglik
  #############################################################################
  if(type == "indep"){
    mc_data_frame <- tidyr::tibble("out" = as.vector(y),
                                   "clust" = as.factor(rep(1:n1, each = n2)))
    form <- out ~ . -clust + (1|clust) + (0 + V1|clust)
  } else if(type == "corr"){
    mc_data_frame <- tidyr::tibble("out" = as.vector(y),
                                   "clust" = as.factor(rep(1:n1, each = n2)))
    form <- out ~ . - clust + (1 + V1|clust)
  } else{
    mc_data_frame <- tidyr::tibble("out" = as.vector(y),
                                   "clust1" = as.factor(rep(1:n1,
                                                            each = n2)),
                                   "clust2" = as.factor(rep(1:n2, n1)))
    form <- out ~ . - clust1 - clust2 + (1|clust1) + (1|clust2)

  }
  # Fit model
  mc_data_frame <- dplyr::bind_cols(mc_data_frame,
                                    tidyr::as_tibble(X[, -1, drop = F]))
  fit <- lme4::lmer(form, data = mc_data_frame, REML = REML)


  # MLEs or REMLEs
  psi_hat <- limestest:::get_psi_hat_lmer(fit)
  Psi_hat <- limestest:::Psi_from_H_cpp(psi_mr = psi_hat[-r], H = H)

  # Get observed information (for Wald) and loglik (for LRT) at estimates
  e <- as.vector(y - X %*% fixef(fit)) # Replace residuals
  Zte <- as.vector(crossprod(Z, e))
  Psi_hat_r <- (1 / psi_hat[r]) * Psi_hat
  stuff_at_est <- limestest:::loglik_psi_cpp(e = e, Z = Z, Zte = Zte, XtZ = XtZ,
                                             ZtZ = ZtZ,
                                             Psi_r = Psi_hat_r,
                                             psi_r = psi_hat[r], H = H,
                                             get_val = TRUE, get_score = FALSE,
                                             get_inf = TRUE, expected = FALSE)
  est_err<- psi_hat - psi
  test_stat_wald <- as.vector(crossprod(est_err, stuff_at_est$inf_mat %*% est_err))

  #############################################################################
  # LR statistics
  #############################################################################
  if(REML){
    ll_null <- stuff_REML$value
    ll_max <-limestest:::res_ll_cpp(Y = y, X = X, Z = Z, XtY = XtY, ZtY = ZtY,
                                         XtX = XtX, XtZ = XtZ, ZtZ = ZtZ,
                                         Psi_r = Psi_hat_r, psi_r = psi_hat[r],
                                         H = H, get_val = TRUE, get_score = FALSE,
                                         get_inf = FALSE)$value
  } else{
    ll_null <- stuff$value
    ll_max  <- stuff_at_est$value
  }

  test_stat_lrt <- 2 * (ll_max - ll_null)

  # Return
  c(test_stat, test_stat_REML, test_stat_lrt, test_stat_wald)
}

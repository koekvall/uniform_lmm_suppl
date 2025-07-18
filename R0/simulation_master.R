library(Matrix)
library(doParallel)
library(lme4)
library(limestest)
library(doRNG)

set.seed(1336)

today <- as.numeric(format(Sys.time(), "%H%d%m%y"))
array_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

###############################################################################
# Settings
###############################################################################
REML <- TRUE # Use REML estimates for competing methods or not
large_n <- FALSE
large_p <- FALSE
# Common settings
num_cores <-  as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
num_reps <- 1e4
num_par <- 10
out_dir <- "/" # has to end in /
fun_dir <- "/" # has to end in /

cl <- makeCluster(num_cores)
registerDoParallel(cl)

source(paste0(fun_dir, "sim_funs.R"))

if(array_id <= num_par){
  jj <- array_id
  setting <- "corr"
} else if(array_id <= 2 * num_par){
  jj <- array_id - num_par
  setting <- "indep"
} else{
  jj <- array_id - 2 * num_par
  setting <- "cross"
}



###############################################################################
# Independent cluster simulation with correlation
###############################################################################

if(setting == "corr"){
  if(large_n){
    n1 <- 1000
    n2 <- 3
    correlations <- seq(-1, -0.8, length.out = num_par)
  } else if(!large_p){
    n1 <- 20
    n2 <- 3
    correlations <- seq(-1, -0.2, length.out = num_par)
  } else{
    n1 <- 20
    n2 <- 3
    correlations <- seq(-1, -0.2, length.out = num_par)
  }

  n <- n1 * n2
  p <- ifelse(large_p, floor(n/5), 2)

  if(large_n & large_p){
      n1 <- 200
      n2 <- 3
      n <- n1 * n2
      p <- 100
      correlations <- seq(-1, -0.4, length.out = num_par)
  }

  psi0 <- 1
  X <- cbind(1, matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p))
  XtX <- crossprod(X)
  b <- rnorm(ncol(X))
  Xb <- X %*% b

  Z <- Matrix::bdiag(lapply(1:n1,
        function(ii)X[((ii - 1) * n2 + 1):(ii * n2), 1:2]))
  ZtZ <- crossprod(Z)
  XtZ <- crossprod(X, Z)
  H1 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(1, 0, 0, 0), 2, 2))
  H2 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 1, 1, 0), 2, 2))
  H3 <- Matrix::kronecker(Matrix::Diagonal(n1), matrix(c(0, 0, 0, 1), 2, 2))
  H <- cbind(H1, H2, H3)
  Psi1 <- matrix(c(1, correlations[jj], correlations[jj], 1), 2, 2)
  Psi1_eig <- eigen(Psi1)
  R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
  R <-  Matrix::kronecker(Matrix::Diagonal(n1), R1)

  out_mat <- foreach(kk = 1:num_reps, .combine = rbind,
                        .errorhandling = "remove",
                        .packages = c("limestest", "lme4", "Matrix")) %dorng%{
                          c(do_one_clust_sim(inner_seed = kk +as.integer(today),
                            n1 = n1,
                            n2 = n2,
                            X = X,
                            Z = Z,
                            Psi1 = Psi1,
                            psi0 = psi0,
                            H = H,
                            R = R,
                            XtX = XtX,
                            XtZ = XtZ,
                            ZtZ = ZtZ,
                            Xb = Xb,
                            REML = REML),
                          n1, n2, p, psi0, correlations[jj], as.numeric(REML))
                        }
}
###############################################################################


###############################################################################
# Independent cluster simulation without correlation
###############################################################################
if(setting == "indep"){
  if(large_n){
    n1 <- 1000
    n2 <- 3
    variances <- seq(0, 0.3, length.out = num_par)
  } else{
    n1 <- 20
    n2 <- 3
    variances <- seq(0, 0.8, length.out = num_par)
  }
  n <- n1 * n2
  p <- ifelse(large_p, floor(sqrt(n)), 2)

  if(large_n & large_p){
      n1 <- 200
      n2 <- 3
      n <- n1 * n2
      p <- 100
      variances <- seq(0, 0.6, length.out = num_par)
  }

  psi0 <- 1
  X <- cbind(1, matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p))
  XtX <- crossprod(X)
  b <- rnorm(ncol(X))
  Xb <- X %*% b
  Z1 <- Matrix::kronecker(Matrix::diag(n1), matrix(1, nrow = n2, ncol = 1))
  Z2 <- Matrix::bdiag(lapply(1:n1,
       function(ii)X[((ii - 1) * n2 + 1):(ii * n2), 2, drop = FALSE]))
  Z <- cbind(Z1, Z2)
  ZtZ <- crossprod(Z)
  XtZ <- crossprod(X, Z)
  H1 <- Matrix::Diagonal(2 * n1, x = c(rep(1, n1), rep(0, n1)))
  H2 <-  Matrix::Diagonal(2 * n1, x = c(rep(0, n1), rep(1, n1)))
  H <- cbind(H1, H2)



  out_mat <- foreach(kk = 1:num_reps, .combine = rbind,
                      .errorhandling = "remove",
                      .packages = c("limestest", "lme4", "Matrix")) %dorng%{
    c(do_one_indep_sim(inner_seed = kk + as.integer(today),
                       n1 = n1,
                       n2 = n2,
                       X = X,
                       Z = Z,
                       Psi = Psi,
                       psi0 = psi0,
                       H = H,
                       R = R,
                       XtX = XtX,
                       XtZ = XtZ,
                       ZtZ = ZtZ,
                       Xb = Xb,
                       REML = REML),
      n1, n2, p, psi0, variances[jj], as.numeric(REML))
  }
}
###############################################################################


###############################################################################
# Crossed random effects simulation
###############################################################################
if(setting == "cross"){
  if(large_n){
    n1 <- 40
    n2 <- 40
    variances <- seq(0, 0.04, length.out = num_par)
  } else{
    n1 <- 10
    n2 <- 10
    variances <- seq(0, 0.5, length.out = num_par) 
  }
  n <- n1 * n2
  p <- ifelse(large_p, floor(n/5), 2)
  
  if(large_n & large_p){
      n1 <- 20
      n2 <- 20
      n <- n1 * n2
      p <- floor(n/5)
      variances <- seq(0, 0.15, length.out = num_par)
  }
  psi0 <- 1
  X <- cbind(1, matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p))
  XtX <- crossprod(X)
  b <- rnorm(ncol(X))
  Xb <- X %*% b
  Z <- cbind(Matrix::kronecker(Matrix::diag(1, nrow = n1), matrix(1, n2, 1)),
                Matrix::kronecker(matrix(1, n1, 1), Matrix::diag(1, nrow = n2)))
  ZtZ <- crossprod(Z)
  XtZ <- crossprod(X, Z)
  H1 <- Matrix::Diagonal(n1 + n2, c(rep(1, n1), rep(0, n2)))
  H2 <- Matrix::Diagonal(n1 + n2, c(rep(0, n1), rep(1, n2)))
  H <- cbind(H1, H2)
  psi <- c(variances[jj], variances[jj])
  Psi <- Matrix::Diagonal(n1 + n2, c(rep(psi[1], n1), rep(psi[2], n2)))

  out_mat <- foreach(kk = 1:num_reps, .combine = rbind,
                      .errorhandling = "remove",
                      .packages = c("limestest", "lme4", "Matrix")) %dorng%{
              c(do_one_cross_sim(inner_seed = kk + as.integer(today),
                n1 = n1,
                n2 = n2,
                X = X,
                Z = Z,
                psi = psi,
                psi0 = psi0,
                H = H,
                XtX = XtX,
                XtZ = XtZ,
                ZtZ = ZtZ,
                Xb = Xb,
                REML = REML),
              n1, n2, p, psi0, variances[jj], as.numeric(REML))
  }
}

###############################################################################

stopCluster(cl)

###############################################################################
# Output
###############################################################################
out <- tidyr::tibble(tidyr::as_tibble(out_mat), "type" = rep(setting, num_reps))
colnames(out) <- c("LL", "RLL", "LRT", "WLD", "seed", "n1", "n2", "p", "psi0",
                  "param", "REML", "type")
saveRDS(out, paste0(out_dir, "lmm_sims_", today, "_", array_id, ".Rds"))

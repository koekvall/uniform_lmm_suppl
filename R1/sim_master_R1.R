library(Matrix)
library(doParallel)
library(lme4)
library(limestest)
library(doRNG)

set.seed(1336)

today <- as.numeric(format(Sys.time(), "%H%d%m%y"))
array_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
num_arrays <- 30 # This needs to be the same as the --array=1-30 argument in sh

###############################################################################
# Settings
###############################################################################

# Common settings
num_cores <-  as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
out_dir <- "/blue/k.ekvall/k.ekvall/Simulations/unif_lmm/R1/Results/" # has to end in /
fun_dir <- "/blue/k.ekvall/k.ekvall/Simulations/unif_lmm/R1/" # has to end in /

cl <- makeCluster(num_cores)
registerDoParallel(cl)

source(paste0(fun_dir, "sim_funs_R1.R"))
# Note the p does not include intercept here; notation slightly different from paper
# to facilitate coding

# psi1 is correlation and psi2 is common variance for correlated random
# intercept and slope
corr_settings <- rbind(cbind("n1" = 1000, "n2" = 3, "p" = 2, "psi1" = 1,
                            "psi2" = seq(-1, -0.8, length.out = 10)),
                      cbind("n1" = 20, "n2" = 3, "p" = 2, "psi1" = 1,
                            "psi2" = seq(-1, -0.2, length.out = 10)),
                      cbind("n1" = 200, "n2" = 3, "p" = 100, "psi1" = 1,
                            "psi2" = seq(-1, -0.4, length.out = 10)))

# psi1 and psi2 are variances of independent random intercept and slope
indep_settings <- rbind(cbind("n1" = 1000, "n2" = 3, "p" = 2,
                             "psi1" = seq(0, 0.3, length.out = 10),
                             "psi2" = seq(0, 0.3, length.out = 10)),
                       cbind("n1" = 1000, "n2" = 3, "p" = 2,
                             "psi1" = seq(0, 0.3, length.out = 10),
                             "psi2" = 0),
                       cbind("n1" = 1000, "n2" = 3, "p" = 2,
                             "psi1" = seq(0, 0.3, length.out = 10),
                             "psi2" = 0.3),
                       cbind("n1" = 20, "n2" = 3, "p" = 2, 
                             "psi1" = seq(0, 0.8, length.out = 10),
                             "psi2" = seq(0, 0.8, length.out = 10)),
                       cbind("n1" = 200, "n2" = 3, "p" = 100,
                             "psi1" = seq(0, 0.6, length.out = 10),
                             "psi2" = seq(0, 0.6, length.out = 10)))

# psi1 and psi2 are variances of independent crossed random effects
cross_settings <- rbind(cbind("n1" = 40, "n2" = 40, "p" = 2,
                             "psi1" = seq(0, 0.04, length.out = 10),
                             "psi2" = seq(0, 0.04, length.out = 10)),
                       cbind("n1" = 20, "n2" = 80, "p" = 2,
                             "psi1" = seq(0, 0.15, length.out = 10),
                             "psi2" = seq(0, 0.15, length.out = 10)),
                       cbind("n1" = 10, "n2" = 10, "p" = 2,
                             "psi1" = seq(0, 0.5, length.out = 10),
                             "psi2" = seq(0, 0.5, length.out = 10)),
                       cbind("n1" = 20, "n2" = 20, "p" = 80,
                             "psi1" = seq(0, 0.15, length.out = 10),
                             "psi2" = seq(0, 0.15, length.out = 10)))
all_settings <- as.data.frame(rbind(corr_settings, indep_settings, cross_settings))
num_corr_set <- nrow(corr_settings)
num_indep_set <- nrow(indep_settings)
num_cross_set <- nrow(cross_settings)
all_settings$type <-  c(rep("corr", num_corr_set), rep("indep", num_indep_set),
                        rep("cross",  num_cross_set))

# Uncomment and edit to run only some of the settings
# all_settings <- all_settings[all_settings$type == "corr", ]

expand_settings <- function(set, seed = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  n <- set$n1 * set$n2
  X <- cbind(1, matrix(runif(n * set$p, min = -1, max = 1), nrow = n, ncol = set$p))
  XtX <- crossprod(X)
  b <- rnorm(set$p)
  Xb <- X %*% b
  
  if(set$type == "corr"){
    psi <- c(set$psi2, set$psi1, set$psi1, set$psi2)
    H1 <- Matrix::kronecker(Matrix::Diagonal(set$n1), matrix(c(1, 0, 0, 0), 2, 2))
    H2 <- Matrix::kronecker(Matrix::Diagonal(set$n1), matrix(c(0, 1, 1, 0), 2, 2))
    H3 <- Matrix::kronecker(Matrix::Diagonal(set$n1), matrix(c(0, 0, 0, 1), 2, 2))
    H <- cbind(H1, H2, H3)
  } else if(set$type == "indep"){
    psi <- c(set$psi1, set$psi2)
    H1 <- Matrix::Diagonal(2 * set$n1, x = c(rep(1, set$n1), rep(0, set$n1)))
    H2 <-  Matrix::Diagonal(2 * set$n1, x = c(rep(0, set$n1), rep(1, set$n1)))
    H <- cbind(H1, H2)
    Z1 <- Matrix::kronecker(Matrix::diag(n1), matrix(1, nrow = n2, ncol = 1))
    Z2 <- Matrix::bdiag(lapply(1:n1,
                               function(ii)X[((ii - 1) * n2 + 1):(ii * n2), 2, drop = FALSE]))
    Z <- cbind(Z1, Z2)
  } else{
    psi <- c(set$psi1, set$psi2)
    H1 <- Matrix::Diagonal(set$n1 + set$n2, c(rep(1, set$n1), rep(0, set$n2)))
    H2 <- Matrix::Diagonal(set$n1 + set$n2, c(rep(0, set$n1), rep(1, set$n2)))
    H <- cbind(H1, H2)
  }
  
  
  list("n1" = set$n1, "n2" = set$n2, "X" = X, "Z" = Z, "psi" = psi, "H" = H, "Xb" = Xb, 
       "R" = NULL, "XtX" = XtX, "XtZ" = XtZ, "ZtZ" = ZtZ, "REML" = TRUE, "type" = set$type)
}


set_per_array <- ceiling(nrow(all_settings)/ num_arrays)
assigned_array <- rep(1:num_arrays, each = set_per_array)[1:nrow(all_settings)]

for(ii in 1:num_run_set){
  if(assigned_array[ii] == array_id){
    
  }
}

  
  psi <- 

  
  XtX <- crossprod(X)
  b <- rnorm(ncol(X))
  Xb <- X %*% b

  Z <- Matrix::bdiag(lapply(1:n1,
        function(ii)X[((ii - 1) * n2 + 1):(ii * n2), 1:2]))
  ZtZ <- crossprod(Z)
  XtZ <- crossprod(X, Z)

  Psi1 <- matrix(c(psi[1:2], psi[2], psi[3]), 2, 2)
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

  psi <-
  X <- cbind(1, matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p))
  XtX <- crossprod(X)
  b <- rnorm(ncol(X))
  Xb <- X %*% b

  ZtZ <- crossprod(Z)
  XtZ <- crossprod(X, Z)

  Psi <- Matrix::Diagonal(2 * n1, x = c(rep(variances[jj], n1),
                                        rep(variances[jj], n1)))
  R <- sqrt(Psi)


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

  Z <- cbind(Matrix::kronecker(Matrix::diag(1, nrow = n1), matrix(1, n2, 1)),
                Matrix::kronecker(matrix(1, n1, 1), Matrix::diag(1, nrow = n2)))
  ZtZ <- crossprod(Z)
  XtZ <- crossprod(X, Z)

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

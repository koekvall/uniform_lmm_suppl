library(Matrix)
library(doParallel)
library(lme4)
library(limestest)
library(doRNG)

set.seed(1336)

today <- as.numeric(format(Sys.time(), "%H%d%m%y"))


###############################################################################
# Settings
###############################################################################

# Common settings
num_reps <- 1e1

# Uncomment this to run on cluster
# num_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
# out_dir <- "/blue/k.ekvall/k.ekvall/Simulations/unif_lmm/R1/Results/" # has to end in /
# fun_dir <- "/blue/k.ekvall/k.ekvall/Simulations/unif_lmm/R1/" # has to end in /
# array_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# num_arrays <- 30 # This needs to be the same as the --array=1-30 argument in sh

# Uncomment this to run on one computer
num_cores <- 8
out_dir <- "~/GitHub/uniform_lmm_suppl/Results/" # has to end in /
fun_dir <- "~/GitHub/uniform_lmm_suppl/" # has to end in /
array_id <- 1
num_arrays <- 1

cl <- makeCluster(num_cores)
registerDoParallel(cl)

source(paste0(fun_dir, "sim_funs_R1.R"))

# psi1-psi3 parameterize 2x2 covariance matrix of random intercept & slope
corr_settings <- rbind(cbind("n1" = 1000, "n2" = 3, "p" = 2, "psi1" = 1,
                            "psi2" = seq(-1, -0.8, length.out = 10),
                            "psi3" = 1, "psi4" = 1),
                      cbind("n1" = 20, "n2" = 3, "p" = 2, "psi1" = 1,
                            "psi2" = seq(-1, -0.2, length.out = 10),
                            "psi3" = 1, "psi4" = 1),
                      cbind("n1" = 200, "n2" = 3, "p" = 100, "psi1" = 1,
                            "psi2" = seq(-1, -0.4, length.out = 10),
                            "psi3" = 1, "psi4" = 1))

# psi1 and psi2 are variances of independent random intercept and slope
indep_settings <- rbind(cbind("n1" = 1000, "n2" = 3, "p" = 2,
                             "psi1" = seq(0, 0.3, length.out = 10),
                             "psi2" = seq(0, 0.3, length.out = 10),
                             "psi3" = 1, "psi4" = NA),
                       cbind("n1" = 1000, "n2" = 3, "p" = 2,
                             "psi1" = seq(0.0001, 0.3, length.out = 10),
                             "psi2" = 0, "psi3" = 1, "psi4" = NA),
                       cbind("n1" = 1000, "n2" = 3, "p" = 2,
                             "psi1" = seq(0, 0.3 - 0.0001, length.out = 10),
                             "psi2" = 0.3, "psi3" = 1, "psi4" = NA),
                       cbind("n1" = 20, "n2" = 3, "p" = 2,
                             "psi1" = seq(0, 0.8, length.out = 10),
                             "psi2" = seq(0, 0.8, length.out = 10),
                             "psi3" = 1, "psi4" = NA),
                       cbind("n1" = 200, "n2" = 3, "p" = 100,
                             "psi1" = seq(0, 0.6, length.out = 10),
                             "psi2" = seq(0, 0.6, length.out = 10),
                             "psi3" = 1, "psi4" = NA))

# psi1 and psi2 are variances of independent crossed random effects
cross_settings <- rbind(cbind("n1" = 40, "n2" = 40, "p" = 2,
                             "psi1" = seq(0, 0.04, length.out = 10),
                             "psi2" = seq(0, 0.04, length.out = 10),
                             "psi3" = 1, "psi4" = NA),
                       cbind("n1" = 20, "n2" = 80, "p" = 2,
                             "psi1" = seq(0, 0.15, length.out = 10),
                             "psi2" = seq(0, 0.15, length.out = 10),
                             "psi3" = 1, "psi4" = NA),
                       cbind("n1" = 10, "n2" = 10, "p" = 2,
                             "psi1" = seq(0, 0.5, length.out = 10),
                             "psi2" = seq(0, 0.5, length.out = 10),
                             "psi3" = 1, "psi4" = NA),
                       cbind("n1" = 20, "n2" = 20, "p" = 80,
                             "psi1" = seq(0, 0.15, length.out = 10),
                             "psi2" = seq(0, 0.15, length.out = 10),
                             "psi3" = 1, "psi4" = NA))
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
  pm1 <- set$p - 1
  X <- cbind(1, matrix(runif(n * pm1, min = -1, max = 1), nrow = n, ncol = pm1))
  XtX <- crossprod(X)
  b <- rnorm(set$p)
  Xb <- X %*% b

  if(set$type == "corr"){
    psi <- c(set$psi1, set$psi2, set$psi3, set$psi4)
    Z <- Matrix::bdiag(lapply(1:set$n1,
                          function(ii)X[((ii - 1) * set$n2 + 1):(ii * set$n2), 1:2]))
    H1 <- Matrix::kronecker(Matrix::Diagonal(set$n1), matrix(c(1, 0, 0, 0), 2, 2))
    H2 <- Matrix::kronecker(Matrix::Diagonal(set$n1), matrix(c(0, 1, 1, 0), 2, 2))
    H3 <- Matrix::kronecker(Matrix::Diagonal(set$n1), matrix(c(0, 0, 0, 1), 2, 2))
    H <- cbind(H1, H2, H3)
    Psi1 <- matrix(psi[c(1, 2, 2, 3)], 2, 2)
    Psi1_eig <- eigen(Psi1)
    R1 <- tcrossprod(diag(sqrt(Psi1_eig$values), ncol(Psi1)), Psi1_eig$vectors)
    R <-  Matrix::kronecker(Matrix::Diagonal(set$n1), R1)
  } else if(set$type == "indep"){
    psi <- c(set$psi1, set$psi2, set$psi3)
    H1 <- Matrix::Diagonal(2 * set$n1, x = c(rep(1, set$n1), rep(0, set$n1)))
    H2 <-  Matrix::Diagonal(2 * set$n1, x = c(rep(0, set$n1), rep(1, set$n1)))
    H <- cbind(H1, H2)
    Z1 <- Matrix::kronecker(Matrix::diag(set$n1),
                            matrix(1, nrow = set$n2, ncol = 1))
    Z2 <- Matrix::bdiag(lapply(1:set$n1,
            function(ii)X[((ii - 1) * set$n2 + 1):(ii * set$n2), 2, drop = FALSE]))
    Z <- cbind(Z1, Z2)
    R <- Matrix::Diagonal(2 * set$n1, x = c(rep(sqrt(psi[1]), set$n1),
                                          rep(sqrt(psi[2]), set$n1)))
  } else{
    psi <- c(set$psi1, set$psi2, set$psi3)
    H1 <- Matrix::Diagonal(set$n1 + set$n2, c(rep(1, set$n1), rep(0, set$n2)))
    H2 <- Matrix::Diagonal(set$n1 + set$n2, c(rep(0, set$n1), rep(1, set$n2)))
    H <- cbind(H1, H2)
    Z <- cbind(Matrix::kronecker(Matrix::Diagonal(n = set$n1), matrix(1, set$n2, 1)),
               Matrix::kronecker(matrix(1, set$n1, 1), Matrix::Diagonal(n = set$n2)))
    R <- methods::as(Matrix::Diagonal(set$n1 + set$n2, c(rep(sqrt(psi[1]), set$n1),
                                             rep(sqrt(psi[2]), set$n2))),
                     "generalMatrix")
  }

  ZtZ <- methods::as(crossprod(Z), "generalMatrix")
  XtZ <- as.matrix(crossprod(X, Z))

  list("n1" = set$n1, "n2" = set$n2, "X" = X, "Z" = Z, "psi" = psi, "H" = H,
       "Xb" = Xb, "R" = R, "XtX" = XtX, "XtZ" = XtZ, "ZtZ" = ZtZ,
       "REML" = 1, "type" = set$type)
}

set_per_array <- ceiling(nrow(all_settings) / num_arrays)
assigned_array <- rep(1:num_arrays, each = set_per_array)[1:nrow(all_settings)]
out <- data.frame()
for(ii in 1:nrow(all_settings)){
  if(assigned_array[ii] == array_id){
    settings_ii <- expand_settings(all_settings[ii, ])
    out_one <- foreach(kk = 1:num_reps, .combine = rbind,
                       .errorhandling = "remove",
                       .packages = c("limestest", "lme4", "Matrix")) %dorng%{
                         merge(data.frame(matrix(do_one_sim(settings_ii), nrow = 1)),
                               all_settings[ii, ])
                       }
    out <- rbind(out, out_one)
  }
}
stopCluster(cl)

###############################################################################
# Output
###############################################################################
names(out)[1:4] <- c("PSCR", "RSCR", "LRT", "WLD")
out <- tidyr::tibble(out_mat)

saveRDS(out, paste0(out_dir, "unif_lmm_sims_", today, "_", array_id, ".Rds"))



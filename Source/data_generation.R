library(tidyverse)
library(SeqNet)
library(ggplot2)
library(MASS)   # For multivariate normal generation
library(stats)  # For Poisson quantile function
library(conquer)
library(quantreg)
library(matrixStats)



## size is the number of subjects, n
## n_vars the number of covariates, p
## for the mean_lambda of poisson
## for the correlation size

generate_multivariate_poisson <- function(n_vars, mean_lambda, rho, size) {
  lambda_shared <- rho * mean_lambda
  lambda_indep <- (1 - rho) * mean_lambda
  shared_poisson <- matrix(rpois(size, lambda_shared), nrow=size, ncol=1)
  indep_poisson <- matrix(rpois(size * n_vars, lambda_indep), nrow=size, ncol=n_vars)
  mv_poisson <- matrix(rep(shared_poisson, n_vars), nrow=size, ncol=n_vars) + indep_poisson
  
  return(mv_poisson)
}



# px is for the number of covariate in x
# mean_lambda: see more details in the multivariate poisson 
# n is for the sample size 
# rho_values: see more details in the multivariate poisson 
# pz is for the number of covariate in z
# nrep is for the number of replication






data_generation <- function(px, mean_lambda, n, rho_values, pz, nrep, seednum) {
  set.seed(seednum)
  for (i in 1:length(n)) {
    for (j in 1:length(rho_values)) {
      for (k in 1:length(pz)) {
        for (rep in 1:nrep) {
          X = generate_multivariate_poisson(px, mean_lambda, rho_values[j], n[i]) 
          Z = matrix(rnorm(n[i] * pz[k]), nrow = n[i], ncol = pz[k])
          Y = matrix(rnorm(n[i]), nrow = n[i])
          data = cbind(Y, Z, X)
          #name = paste("simu_data/n_", n[i], "_rho_", rho_values[j], "_pz_", pz[k], "_rep_", rep, ".csv", sep = "")
        }
      }
    }
  }
}


data_generation2 <- function(px, mean_lambda, n, rho_values, pz, nrep, seednum) {
  set.seed(seednum)
  p <- runif(px, 0, 0.5)
  for (i in 1:length(n)) {
    for (k in 1:length(pz)) {
      for (rep in 1:nrep) {
        X = sapply(p, function(pi) rbinom(N, size = 2, prob = pi))
        X = matrix(X, nrow = n[i], ncol = px)
        Z = matrix(rnorm(n[i] * pz[k]), nrow = n[i], ncol = pz[k])
        Y = matrix(rnorm(n[i]), nrow = n[i])
        data = cbind(Y, Z, X)
      }
    }
  }
}




data_generation3 <- function(n, pz, nrep) {
  for (i in 1:length(n)) {
    for (k in 1:length(pz)) {
      for (rep in 1:nrep) {
        X1 = rbinom(n[i], 2, 0.5)
        X2 = rbinom(n[i], 2, 0.5)
        rho <- 0.3  # Correlation between Z variables
        Sigma_Z <- matrix(rho, nrow = pz[k], ncol = pz[k])  # Fill with rho
        diag(Sigma_Z) <- 1
        Z = mvrnorm(n = n[i], mu = rep(0, pz[k]), Sigma = Sigma_Z)
        error = matrix(rt(n[i], df = 2), nrow = n[i])
        Y = 0.5 + Z %*% matrix(rep(0.002, pz[k]), ncol = 1) + 0.3 * X1 + 0.4 * X2 + (3 + 0.5 * X1 + 0.5 * X2) * error
        data = cbind(Y, Z, X1, X2)
      }
    }
  }
  return(data)
}



data_generation4 <- function(n = 300, p = p, z_dim = 30, x_block_size = x_block_size, 
                              num_blocks = num_blocks, maf = 0.3, k = 1, seednum) {
  set.seed(seednum)
  omega_k = 0.1 * k
  Z = matrix(rnorm(n * z_dim), nrow = n, ncol = z_dim)
  X = matrix(rbinom(n * p, size = 2, prob = maf), nrow = n, ncol = p)
  Y_list = vector("list", num_blocks)
  for (b in 1:1) {
    x_idx = ((b - 1) * x_block_size + 1):(b * x_block_size)
    X_block = X[, x_idx]
    if (x_block_size == 1) {
      signal_mean = 0.5 + Z %*% rep(0.2, z_dim) + omega_k * sum(X_block)
      hetero_sd = 3 + X_block * rep(0.2, x_block_size)
    } else {
      X_block = X[, x_idx]
      signal_mean = 0.5 + Z %*% rep(0.2, z_dim) + omega_k * rowSums(X_block)
      hetero_sd = 3 + X_block %*% rep(0.2, x_block_size)
    }
    epsilon = rnorm(n)
    Y = signal_mean + hetero_sd * epsilon
  }
  data = cbind(Y, Z, X)
  return(data)
}


bias <- function(n, replications, tau, error_distribution, p = 10, seednum = 123) {
  set.seed(seednum)
  beta = c(rep(0.5, p / 2), rep(0.7, p / 2))
  gamma = c(rep(c(-0.1, 0.1), p / 2))
  if (error_distribution == "normal") {
    z_tau = qnorm(tau)
  } else if (error_distribution == "cauchy") {
    z_tau = qcauchy(tau)
  }
  beta_tau_true = beta + z_tau * gamma
  
  bias_qr_vec = numeric(replications)
  bias_conquer_vec = numeric(replications)
  
  
  for (r in 1:replications) {
    X = matrix(rnorm(n * p), n, p)
    scale = 0.3 + X %*% gamma
    eps = rnorm(n)
    y = 1 + X %*% beta + scale * eps
    
    qr_fit = rq(y ~ X, tau = tau)
    conquer_fit = conquer(X, y, tau = tau, kernel = "Gaussian")
    
    coef_qr = coef(qr_fit)
    coef_conquer = conquer_fit$coeff
    bias_qr_vec[r] = mean(coef_qr - beta_tau_true)
    bias_conquer_vec[r] = mean(coef_conquer - beta_tau_true)
  }
  bias_qr = mean(bias_qr_vec)
  bias_conquer = mean(bias_conquer_vec)
  return(c(bias_qr, bias_conquer))
}

# Pseudo Marginal NPP Algorithm
# Core functions for pseudo-marginal MCMC with normalized power priors

library(MASS)
library(mvtnorm)
library(foreach)
library(doParallel)

# -----------------------------
# Model utilities
# -----------------------------

sample_theta_power <- function(n, tau, dbar, invSigma, invSigma0, N0) {
  Q <- tau * N0 * invSigma + invSigma0
  V <- solve(Q)
  mu <- V %*% (tau * N0 * invSigma %*% dbar)
  mvrnorm(n, mu, V)
}

log_lik0 <- function(D0, theta, Sigma) {
  sum(dmvnorm(D0, mean = theta, sigma = Sigma, log = TRUE))
}

log_prior_theta <- function(theta, Sigma0) {
  p <- length(theta)
  sum(dmvnorm(theta, mean = rep(0, p), sigma = Sigma0, log = TRUE))
}

log_lik_data <- function(D, theta, Sigma) {
  sum(dmvnorm(D, mean = theta, sigma = Sigma, log = TRUE))
}

log_prior_eta <- function(eta, alpha0, beta0) {
  dbeta(eta, alpha0, beta0, log = TRUE)
}

propose_theta <- function(theta, sigma_prop = 0.5) {
  p <- length(theta)
  as.numeric(mvrnorm(1, theta, sigma_prop^2 * diag(p)))
}

propose_eta <- function(eta, sigma_prop = 0.05) {
  x <- runif(1, eta - sigma_prop, eta + sigma_prop)
  if (x < 0 || x > 1) eta else x
}

# -----------------------------
# Exact reference (debug only)
# -----------------------------

log_c0 <- function(eta, dbar, Sigma, Sigma0, N0, S) {
  p <- length(dbar)
  Sigma_inv <- solve(Sigma)
  Sigma0_inv <- solve(Sigma0)

  A <- eta * N0 * Sigma_inv + Sigma0_inv
  A_inv <- solve(A)
  B <- eta * N0 * Sigma_inv %*% dbar

  term1 <- -eta * N0 * p / 2 * log(2 * pi)
  term2 <- -eta * N0 / 2 * log(det(Sigma))
  term3 <- -eta / 2 * sum(diag(S %*% Sigma_inv))
  term4 <- -eta * N0 / 2 * t(dbar) %*% Sigma_inv %*% dbar
  term5 <- t(B) %*% A_inv %*% B / 2
  term6 <- -0.5 * log(det(Sigma0)) - 0.5 * log(det(A))

  as.numeric(term1 + term2 + term3 + term4 + term5 + term6)
}

unbiased_ratio_exact <- function(eta0, eta1, dbar, Sigma, Sigma0, N0, S) {
  if (eta0 == eta1) return(1.0)
  exp(log_c0(eta1, dbar, Sigma, Sigma0, N0, S) -
        log_c0(eta0, dbar, Sigma, Sigma0, N0, S))
}

# -----------------------------
# Single chain kernel
# -----------------------------

run_single_mcmc_chain <- function(chain_id, niter, init_theta, init_eta,
                                  D, D0, Sigma, Sigma0,
                                  alpha0, beta0,
                                  logratio_estimator,
                                  sigma_prop_theta = 0.5,
                                  sigma_prop_eta = 0.05,
                                  use_exact = FALSE,
                                  seed = NULL,
                                  verbose = TRUE) {

  if (!is.null(seed)) set.seed(seed + chain_id)

  p <- length(init_theta)
  N0 <- nrow(D0)

  dbar <- colMeans(D0)
  invSigma <- solve(Sigma)
  invSigma0 <- solve(Sigma0)
  S <- t(D0 - matrix(dbar, N0, p, byrow = TRUE)) %*%
       (D0 - matrix(dbar, N0, p, byrow = TRUE))

  chain_th <- matrix(NA, nrow = niter, ncol = p)
  chain_et <- numeric(niter)
  chain_th[1,] <- init_theta
  chain_et[1] <- init_eta

  accept_theta <- 0
  accept_eta <- 0

  for (i in 2:niter) {

    # theta update
    th0 <- chain_th[i-1,]
    et0 <- chain_et[i-1]
    th1 <- propose_theta(th0, sigma_prop_theta)

    lp0 <- log_prior_theta(th0, Sigma0) +
           log_lik_data(D, th0, Sigma) +
           et0 * log_lik0(D0, th0, Sigma)

    lp1 <- log_prior_theta(th1, Sigma0) +
           log_lik_data(D, th1, Sigma) +
           et0 * log_lik0(D0, th1, Sigma)

    if (log(runif(1)) < lp1 - lp0) {
      chain_th[i,] <- th1
      accept_theta <- accept_theta + 1
    } else {
      chain_th[i,] <- th0
    }

    # eta update
    et0 <- chain_et[i-1]
    th_curr <- chain_th[i,]
    et1 <- propose_eta(et0, sigma_prop_eta)

    lpa0 <- log_prior_eta(et0, alpha0, beta0)
    lpa1 <- log_prior_eta(et1, alpha0, beta0)

    hist_ratio <- (et1 - et0) * log_lik0(D0, th_curr, Sigma)

    if (use_exact) {
      logRhat <- log(unbiased_ratio_exact(et0, et1, dbar, Sigma, Sigma0, N0, S))
    } else {
      logRhat <- logratio_estimator(
        et0, et1,
        D0 = D0,
        Sigma = Sigma,
        dbar = dbar,
        invSigma = invSigma,
        invSigma0 = invSigma0,
        N0 = N0,
        S = S
      )
    }

    log_r <- (lpa1 - lpa0) + hist_ratio - logRhat

    if (log(runif(1)) < log_r) {
      chain_et[i] <- et1
      accept_eta <- accept_eta + 1
    } else {
      chain_et[i] <- et0
    }

    if (verbose && i %% 500 == 0) {
      cat("Chain", chain_id, "iter", i,
          "| acc theta =", round(accept_theta / i, 3),
          "| acc eta =", round(accept_eta / i, 3), "\n")
    }
  }

  list(theta = chain_th,
       eta = chain_et,
       accept_rates = c(theta = accept_theta / niter,
                        eta = accept_eta / niter))
}

# -----------------------------
# Multiple chains
# -----------------------------

run_joint_mcmc <- function(niter, init_theta, init_eta,
                           D, D0, Sigma, Sigma0,
                           alpha0, beta0,
                           logratio_estimator,
                           sigma_prop_theta = 0.5,
                           sigma_prop_eta = 0.05,
                           chains = 4,
                           iter_warmup = 1000,
                           seed = 123,
                           verbose = TRUE,
                           use_exact = FALSE) {

  total_iter <- niter

  init_vals <- vector("list", chains)
  p <- length(init_theta)

  for (i in 1:chains) {
    set.seed(seed + i)
    init_vals[[i]] <- list(
      theta = init_theta + rnorm(p, 0, 0.1),
      eta   = max(0.01, min(0.99, init_eta + rnorm(1, 0, 0.05)))
    )
  }

  chain_results <- vector("list", chains)

  for (chain_id in 1:chains) {
    chain_results[[chain_id]] <- run_single_mcmc_chain(
      chain_id = chain_id,
      niter = total_iter,
      init_theta = init_vals[[chain_id]]$theta,
      init_eta   = init_vals[[chain_id]]$eta,
      D = D, D0 = D0, Sigma = Sigma, Sigma0 = Sigma0,
      alpha0 = alpha0, beta0 = beta0,
      logratio_estimator = logratio_estimator,
      sigma_prop_theta = sigma_prop_theta,
      sigma_prop_eta   = sigma_prop_eta,
      use_exact = use_exact,
      seed = seed,
      verbose = verbose
    )
  }

  list(
    chain_results = chain_results,
    iter_warmup = iter_warmup,
    chains = chains
  )
}

# -----------------------------
# Draws extractor
# -----------------------------

mcmc_to_draws_df <- function(mcmc_results) {

  chain_results <- mcmc_results$chain_results
  iter_warmup <- mcmc_results$iter_warmup
  chains <- length(chain_results)
  p <- ncol(chain_results[[1]]$theta)

  all_samples <- list()

  for (chain_id in 1:chains) {

    theta_samples <- chain_results[[chain_id]]$theta[-(1:iter_warmup), , drop = FALSE]
    eta_samples   <- chain_results[[chain_id]]$eta[-(1:iter_warmup)]

    n_samples <- nrow(theta_samples)

    chain_df <- data.frame(
      .chain = rep(chain_id, n_samples),
      .iteration = 1:n_samples,
      .draw = (chain_id - 1) * n_samples + 1:n_samples,
      eta = eta_samples
    )

    for (j in 1:p) {
      chain_df[[paste0("theta[", j, "]")]] <- theta_samples[, j]
    }

    all_samples[[chain_id]] <- chain_df
  }

  draws_df <- do.call(rbind, all_samples)
  class(draws_df) <- c("draws_df", "draws", "tbl_df", "tbl", "data.frame")
  draws_df
}

# -----------------------------
# High-level wrapper
# -----------------------------

run_pseudo_marginal <- function(scn,
                                logratio_estimator,
                                n_iter = 3000,
                                chains = 2,
                                iter_warmup = 500,
                                sigma_prop_theta = 0.5,
                                sigma_prop_eta = 0.05,
                                seed = 123,
                                verbose = TRUE,
                                use_exact = FALSE) {

  init_theta <- rep(0, scn$p)
  init_eta   <- 0.5

  fit <- run_joint_mcmc(
    niter = n_iter,
    init_theta = init_theta,
    init_eta   = init_eta,
    D = scn$d,
    D0 = scn$d0,
    Sigma = scn$Sigma,
    Sigma0 = scn$Sigma0,
    alpha0 = 1.0,
    beta0  = 1.0,
    logratio_estimator = logratio_estimator,
    sigma_prop_theta = sigma_prop_theta,
    sigma_prop_eta   = sigma_prop_eta,
    chains = chains,
    iter_warmup = iter_warmup,
    seed = seed,
    verbose = verbose,
    use_exact = use_exact
  )

  draws_df <- mcmc_to_draws_df(fit)

  list(
    draws = draws_df,
    eta   = draws_df$eta,
    fit   = fit
  )
}
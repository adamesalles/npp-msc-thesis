# estimators.R
# Wrappers for marginal likelihood estimators

library(mvtnorm)

# log L0 for MVN
loglik0 <- function(theta, dbar, invSigma, N0) {
  quad <- rowSums((theta - dbar) * (invSigma %*% t(theta - dbar)))
  -0.5 * N0 * quad
}

# prior estimator
delta_prior <- function(theta, eta, dbar, invSigma, N0) {
  mean(exp(eta * loglik0(theta, dbar, invSigma, N0)))
}

# Importance sampling around MLE
delta_is <- function(theta, eta, dbar, invSigma, N0, log_w) {
  mean(exp(eta * loglik0(theta, dbar, invSigma, N0) + log_w))
}

# Single-sample log-ratio estimator
log_ratio_single <- function(eta0, eta1, theta, dbar, invSigma, N0) {
  (eta1 - eta0) * loglik0(theta, dbar, invSigma, N0)
}

# Numerically stable log(mean(exp(x))) using the log-sum-exp trick
log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Log-ratio via prior_importance estimator
make_logratio_prior_importance <- function(K = 200, n = 200) {
  function(eta0, eta1, D0, Sigma, dbar, invSigma, invSigma0, N0, S) {
    if (eta0 == eta1) return(0.0)

    theta <- mvtnorm::rmvnorm(K, mean = rep(0, length(dbar)), sigma = solve(invSigma0))
    loglik <- apply(theta, 1, function(th) {
      sum(dmvnorm(D0, mean = th, sigma = Sigma, log = TRUE))
    })

    log_mean_exp(eta1 * loglik) - log_mean_exp(eta0 * loglik)
  }
}

# Log-ratio via importance sampling around MLE
make_logratio_importance <- function(theta_mle, Sigma_mle, K = 200) {
  function(eta0, eta1, D0, Sigma, dbar, invSigma, invSigma0, N0, S) {
    if (eta0 == eta1) return(0.0)

    theta <- mvtnorm::rmvnorm(K, mean = theta_mle, sigma = Sigma_mle)

    loglik <- apply(theta, 1, function(th) {
      sum(dmvnorm(D0, mean = th, sigma = Sigma, log = TRUE))
    })

    logw <- apply(theta, 1, function(th) {
      dmvnorm(th, mean = rep(0, length(th)), sigma = solve(invSigma0), log = TRUE) -
        dmvnorm(th, mean = theta_mle, sigma = Sigma_mle, log = TRUE)
    })

    log_mean_exp(eta1 * loglik + logw) - log_mean_exp(eta0 * loglik + logw)
  }
}

# Log-ratio via power posterior Monte Carlo integration
make_logratio_power <- function(K = 50, n = 50) {
  function(eta0, eta1, D0, Sigma, dbar, invSigma, invSigma0, N0, S) {
    if (eta0 == eta1) return(0.0)

    etl <- min(eta0, eta1)
    etu <- max(eta0, eta1)
    tau <- runif(K, etl, etu)

    est <- numeric(K)

    for (k in seq_len(K)) {
      theta_star <- sample_theta_power(n, tau[k], dbar, invSigma, invSigma0, N0)
      ll <- apply(theta_star, 1, function(th) {
        sum(dmvnorm(D0, mean = th, sigma = Sigma, log = TRUE))
      })
      est[k] <- mean(ll)
    }

    Delta <- mean(est) * (eta1 - eta0)
    Delta
  }
}
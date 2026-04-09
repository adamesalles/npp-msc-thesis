# scenario.R
# Data generation for discrepancy regimes

library(MASS)

make_scenario <- function(p = 3,
                          N0 = 50,
                          N = 50,
                          discrepancy = c("small","medium","large"),
                          seed = 1) {

  set.seed(seed)
  discrepancy <- match.arg(discrepancy)

  Sigma  <- diag(p)
  Sigma0 <- diag(p) * 2

  theta_true <- rep(0, p)

  shift <- switch(
    discrepancy,
    small  = rep(0.1, p),
    medium = rep(0.5, p),
    large  = rep(1.0, p)
  )

  d0 <- MASS::mvrnorm(N0, theta_true + shift, Sigma)
  d  <- MASS::mvrnorm(N,  theta_true, Sigma)

  list(
    p = p,
    N0 = N0,
    N  = N,
    Sigma  = Sigma,
    Sigma0 = Sigma0,
    d0 = d0,
    d  = d,
    dbar0 = colMeans(d0)
  )
}
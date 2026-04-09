# baseline_stan.R
# Baseline Stan fit

library(cmdstanr)
library(posterior)

run_baseline <- function(scn,
                         iter = 10000,
                         chains = 4,
                         alpha0 = 1,
                         beta0 = 1) {

  mod <- cmdstan_model("../npp_mvn.stan")

  dat <- list(
    p = scn$p,
    N0 = scn$N0,
    N  = scn$N,
    D0 = scn$d0,
    D  = scn$d,
    Sigma  = scn$Sigma,
    Sigma0 = scn$Sigma0,
    alpha0 = alpha0,
    beta0  = beta0
  )

  fit <- mod$sample(
    data = dat,
    iter_sampling = iter,
    iter_warmup = iter / 2,
    chains = chains,
    refresh = 0
  )

  draws <- as_draws_df(fit$draws("eta"))

  list(
    fit = fit,
    eta = draws$eta,
    summary = fit$summary("eta")
  )
}
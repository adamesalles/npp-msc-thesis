# run_all.R
# End-to-end experiment runner

library(dplyr)
library(tidyr)
library(purrr)
library(readr)

source("scenario.R")
source("estimators.R")
source("baseline_stan.R")
source("summarize.R")
source("plots.R")
source("pseudo_marginal_algorithm.R")

# -------------------------------
# Quick sanity test
# -------------------------------
RUN_TEST <- FALSE

if (RUN_TEST) {

  cat("Running sanity test...\n")

  scn <- make_scenario(discrepancy = "small", seed = 1)


  # Stan baseline
  base <- run_baseline(scn)
  eta_ref <- base$eta

  # Exact reference
  pm_exact <- run_pseudo_marginal(
    scn,
    logratio_estimator = NULL,
    use_exact = TRUE,
    n_iter = 5000,
    chains = 4,
    iter_warmup = 1000
  )

  # prior_importance
  logratio_prior_importance <- make_logratio_prior_importance(K = 50)
  pm_prior_importance <- run_pseudo_marginal(
    scn,
    logratio_estimator = logratio_prior_importance,
    n_iter = 5000,
    chains = 4,
    iter_warmup = 1000
  )

  # Importance
  theta_mle <- colMeans(scn$d0)
  Sigma_mle <- cov(scn$d0) + 1e-6 * diag(scn$p)

  logratio_is <- make_logratio_importance(
    theta_mle = theta_mle,
    Sigma_mle = Sigma_mle,
    K = 50
  )

  pm_is <- run_pseudo_marginal(
    scn,
    logratio_estimator = logratio_is,
    n_iter = 5000,
    chains = 4,
    iter_warmup = 1000
  )

  # Power posterior
  logratio_power <- make_logratio_power(K = 10, n = 10)

  pm_power <- run_pseudo_marginal(
    scn,
    logratio_estimator = logratio_power,
    n_iter = 5000,
    chains = 4,
    iter_warmup = 1000
  )

  # cat("\nStan baseline summary:\n")
  print(summary(eta_ref))
  print(summary(pm_exact$eta))
  print(summary(pm_prior_importance$eta))
  print(summary(pm_is$eta))
  print(summary(pm_power$eta))

  df <- bind_rows(
    tibble(eta = eta_ref,       method = "stan"),
    tibble(eta = pm_exact$eta,  method = "exact"),
    tibble(eta = pm_prior_importance$eta, method = "prior_importance"),
    tibble(eta = pm_is$eta,     method = "importance"),
    tibble(eta = pm_power$eta,  method = "power")
  )

  p <- plot_posteriors(df, "Sanity test posterior comparison")
  print(p)

  stop("Sanity test finished. Set RUN_TEST <- FALSE for full runs.")
}

# -------------------------------
# Full experiments
# -------------------------------

set.seed(123)
scenarios <- c("small","medium","large")

for (sc in scenarios) {

  cat("Running scenario:", sc, "\n")
  scn <- make_scenario(discrepancy = sc)

  base <- run_baseline(scn)
  eta_ref <- base$eta

  # prior_importance
  logratio_prior_importance <- make_logratio_prior_importance(K = 100)
  pm_prior_importance <- run_pseudo_marginal(
    scn,
    logratio_estimator = logratio_prior_importance,
    n_iter = 5000,
    chains = 4,
    iter_warmup = 1000
  )

  # Power posterior
  logratio_power <- make_logratio_power(K = 10, n = 10)

  pm_power <- run_pseudo_marginal(
    scn,
    logratio_estimator = logratio_power,
    n_iter = 5000,
    chains = 4,
    iter_warmup = 1000
  )

  # Importance
  theta_mle <- colMeans(scn$d0)
  Sigma_mle <- cov(scn$d0) + 1e-6 * diag(scn$p)

  logratio_is <- make_logratio_importance(
    theta_mle = theta_mle,
    Sigma_mle = Sigma_mle,
    K = 100
  )

  pm_is <- run_pseudo_marginal(
    scn,
    logratio_estimator = logratio_is,
    n_iter = 5000,
    chains = 4,
    iter_warmup = 1000
  )

  df_plot <- bind_rows(
    tibble(eta = eta_ref, method = "baseline"),
    tibble(eta = pm_prior_importance$eta, method = "prior_importance"),
    tibble(eta = pm_power$eta, method = "power"),
    tibble(eta = pm_is$eta, method = "importance")
  )

  p <- plot_posteriors(df_plot, paste("Discrepancy:", tools::toTitleCase(sc)))
  ggsave(paste0("posterior_", sc, ".png"), p, width = 7, height = 5)

  tab <- bind_rows(
    summarize_vs_baseline(pm_prior_importance$eta, eta_ref) %>% mutate(method = "prior_importance"),
    summarize_vs_baseline(pm_is$eta,   eta_ref) %>% mutate(method = "importance"),
    summarize_vs_baseline(pm_power$eta,   eta_ref) %>% mutate(method = "power"),
    summarize_vs_baseline(eta_ref,  eta_ref) %>% mutate(method = "baseline")
  )

  readr::write_csv(tab, paste0("summary_", sc, ".csv"))
}
# run_convergence.R
# Convergence-based experiment runner for NPP pseudo-marginal MCMC.
#
# For each scenario (small/medium/large discrepancy), every estimator —
# including the exact one — is run in independent batches of
# chains_per_batch × n_iter_per_batch iterations until its KS distance
# to the Stan baseline drops below ks_threshold.
#
# "Exact" refers to the pseudo-marginal algorithm using the closed-form
# normalizing constant ratio (log_c0), which targets the true posterior.
# It serves as a reference: it should converge fastest, and any estimator
# that needs significantly more chains/time to converge is a weaker choice.
#
# Checkpoints overwrite each batch (no numbered subfolders):
#   {output_dir}/
#     convergence.csv    — cumulative log: batch, total_chains, elapsed_sec,
#                          ks_stat, converged, eta_mean, eta_sd
#     posterior.png      — accumulated samples vs baseline (overwritten)
#     traces/            — trace plots for the most recent batch (overwritten)
#       chain_01_eta.png
#       ...
#
# Final outputs per scenario:
#   results/{sc}/comparison_summary.csv
#   results/{sc}/final_comparison.png

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)

source("scenario.R")
source("estimators.R")
source("baseline_stan.R")
source("summarize.R")
source("plots.R")
source("pseudo_marginal_algorithm.R")

# -----------------------------
# Convergence runner
# -----------------------------

run_until_convergence <- function(scn,
                                   logratio_estimator,
                                   estimator_name,
                                   eta_ref,
                                   output_dir,
                                   use_exact        = FALSE,
                                   ks_threshold     = 0.05,
                                   n_iter_per_batch = 5000,
                                   iter_warmup      = 1000,
                                   chains_per_batch = 4,
                                   max_batches      = 10) {

  dir.create(file.path(output_dir, "traces"),
             recursive = TRUE, showWarnings = FALSE)

  all_eta_samples <- numeric(0)
  convergence_log <- list()
  start_time      <- proc.time()["elapsed"]

  for (batch_num in seq_len(max_batches)) {

    cat(sprintf("[%s] Batch %d / %d ...\n", estimator_name, batch_num, max_batches))

    seed_batch <- 123 + (batch_num - 1) * chains_per_batch

    batch_fit <- run_joint_mcmc(
      niter              = n_iter_per_batch,
      init_theta         = rep(0, scn$p),
      init_eta           = 0.5,
      D                  = scn$d,
      D0                 = scn$d0,
      Sigma              = scn$Sigma,
      Sigma0             = scn$Sigma0,
      alpha0             = 1.0,
      beta0              = 1.0,
      logratio_estimator = logratio_estimator,
      sigma_prop_theta   = 0.5,
      sigma_prop_eta     = 0.05,
      chains             = chains_per_batch,
      iter_warmup        = iter_warmup,
      seed               = seed_batch,
      verbose            = TRUE,
      use_exact          = use_exact
    )

    # Extract post-warmup eta for each chain in this batch
    batch_eta_by_chain <- lapply(
      seq_len(chains_per_batch),
      function(i) batch_fit$chain_results[[i]]$eta[-(1:iter_warmup)]
    )

    all_eta_samples <- c(all_eta_samples, unlist(batch_eta_by_chain))

    total_chains <- batch_num * chains_per_batch
    elapsed_sec  <- as.numeric(proc.time()["elapsed"] - start_time)
    ks_stat      <- as.numeric(suppressWarnings(
      ks.test(all_eta_samples, eta_ref)$statistic
    ))
    converged <- ks_stat < ks_threshold

    convergence_log[[batch_num]] <- list(
      batch        = batch_num,
      total_chains = total_chains,
      elapsed_sec  = elapsed_sec,
      ks_stat      = ks_stat,
      converged    = converged,
      eta_mean     = mean(all_eta_samples),
      eta_sd       = sd(all_eta_samples)
    )

    cat(sprintf(
      "[%s] Batch %d: chains=%d, elapsed=%.1fs, KS=%.4f, mean=%.3f, sd=%.3f%s\n",
      estimator_name, batch_num, total_chains, elapsed_sec, ks_stat,
      mean(all_eta_samples), sd(all_eta_samples),
      if (converged) " -> CONVERGED" else ""
    ))

    # --- convergence.csv (cumulative, overwrites) ---
    conv_df <- do.call(rbind, lapply(convergence_log, as.data.frame))
    readr::write_csv(conv_df, file.path(output_dir, "convergence.csv"))

    # --- posterior.png (overwrites) ---
    legend_label <- sprintf(
      "%s [chains=%d, time=%.1fs, KS=%.3f]",
      estimator_name, total_chains, elapsed_sec, ks_stat
    )
    df_plot <- dplyr::bind_rows(
      tibble::tibble(eta = eta_ref,         method = "baseline"),
      tibble::tibble(eta = all_eta_samples, method = legend_label)
    )
    p_post <- plot_posteriors(df_plot, title = estimator_name)
    ggplot2::ggsave(file.path(output_dir, "posterior.png"), p_post,
                   width = 7, height = 5)

    # --- traces/ (overwrites — shows mixing of the most recent batch) ---
    for (i in seq_len(chains_per_batch)) {
      p_trace <- plot_trace(
        eta_chain = batch_eta_by_chain[[i]],
        chain_id  = i,
        title     = sprintf("%s | batch %d, chain %d",
                            estimator_name, batch_num, i)
      )
      ggplot2::ggsave(
        file.path(output_dir, "traces", sprintf("chain_%02d_eta.png", i)),
        p_trace, width = 7, height = 3
      )
    }

    if (converged) break
  }

  final_log <- convergence_log[[length(convergence_log)]]

  list(
    eta          = all_eta_samples,
    total_chains = final_log$total_chains,
    elapsed_sec  = final_log$elapsed_sec,
    batches_run  = length(convergence_log),
    converged    = final_log$converged,
    final_ks     = final_log$ks_stat
  )
}

# -----------------------------
# Main experiment loop
# -----------------------------

set.seed(123)
scenarios <- c("small", "medium", "large")

for (sc in scenarios) {

  cat(sprintf("\n===== Scenario: %s =====\n", sc))

  results_dir <- file.path("results", sc)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  scn <- make_scenario(discrepancy = sc)

  # Stan baseline (fixed reference — not iterated)
  cat("Running Stan baseline...\n")
  base    <- run_baseline(scn)
  eta_ref <- base$eta

  # Importance estimator needs MLE precomputed from historical data
  theta_mle <- colMeans(scn$d0)
  Sigma_mle <- cov(scn$d0) + 1e-6 * diag(scn$p)

  # All estimators run through the same convergence loop.
  # "exact" uses the closed-form normalizing constant ratio (use_exact = TRUE,
  # logratio_estimator = NULL). It should converge fastest and is the reference
  # for how many chains each approximate estimator needs in comparison.
  estimators_list <- list(
    exact            = list(logratio_estimator = NULL,   use_exact = TRUE),
    prior_importance = list(logratio_estimator = make_logratio_prior_importance(K = 100),
                            use_exact = FALSE),
    importance       = list(logratio_estimator = make_logratio_importance(
                              theta_mle = theta_mle,
                              Sigma_mle = Sigma_mle,
                              K = 100
                            ), use_exact = FALSE),
    power            = list(logratio_estimator = make_logratio_power(K = 10, n = 10),
                            use_exact = FALSE)
  )

  conv_results <- list()

  for (est_name in names(estimators_list)) {
    cat(sprintf("\nRunning estimator: %s\n", est_name))
    est <- estimators_list[[est_name]]
    conv_results[[est_name]] <- run_until_convergence(
      scn                = scn,
      logratio_estimator = est$logratio_estimator,
      estimator_name     = est_name,
      eta_ref            = eta_ref,
      output_dir         = file.path(results_dir, est_name),
      use_exact          = est$use_exact
    )
  }

  # Comparison summary CSV
  comparison_df <- dplyr::bind_rows(lapply(names(conv_results), function(nm) {
    r <- conv_results[[nm]]
    tibble::tibble(
      estimator    = nm,
      total_chains = r$total_chains,
      elapsed_sec  = r$elapsed_sec,
      batches_run  = r$batches_run,
      final_ks     = r$final_ks,
      converged    = r$converged
    )
  }))
  readr::write_csv(comparison_df,
                   file.path(results_dir, "comparison_summary.csv"))
  cat("Saved comparison_summary.csv\n")

  # Final comparison plot: all estimators + baseline
  make_label <- function(nm, r) {
    sprintf("%s [chains=%d, time=%.1fs, KS=%.3f]",
            nm, r$total_chains, r$elapsed_sec, r$final_ks)
  }

  df_final <- dplyr::bind_rows(
    tibble::tibble(eta = eta_ref, method = "baseline"),
    purrr::map_dfr(names(conv_results), function(nm) {
      tibble::tibble(eta    = conv_results[[nm]]$eta,
                     method = make_label(nm, conv_results[[nm]]))
    })
  )

  p_final <- plot_posteriors(
    df_final,
    title = paste("Discrepancy:", tools::toTitleCase(sc), "— final comparison")
  )
  ggplot2::ggsave(file.path(results_dir, "final_comparison.png"),
                 p_final, width = 9, height = 7)
  cat("Saved final_comparison.png\n")
}

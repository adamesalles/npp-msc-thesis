# plots.R
# Posterior plots

library(ggplot2)

plot_posteriors <- function(df, title) {
  ggplot(df, aes(x = eta, color = method)) +
    geom_density(linewidth = 1) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(ncol = 1)) +
    labs(title = title, x = "eta", y = "density")
}

plot_trace <- function(eta_chain, chain_id, title) {
  df <- data.frame(iteration = seq_along(eta_chain), eta = eta_chain)
  ggplot(df, aes(x = iteration, y = eta)) +
    geom_line(linewidth = 0.4, colour = "steelblue") +
    theme_minimal(base_size = 14) +
    labs(title = title, x = "Iteration (post-warmup)", y = "eta")
}
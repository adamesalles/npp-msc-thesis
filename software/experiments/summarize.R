# summarize.R
# Quantitative summaries

library(posterior)

summarize_vs_baseline <- function(x, baseline) {
  ks <- suppressWarnings(ks.test(x, baseline)$statistic)

  tibble::tibble(
    mean = mean(x),
    sd = sd(x),
    q05 = quantile(x, 0.05),
    q50 = quantile(x, 0.50),
    q95 = quantile(x, 0.95),
    ks = as.numeric(ks)
  )
}
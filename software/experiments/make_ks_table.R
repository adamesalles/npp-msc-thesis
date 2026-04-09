# make_ks_table.R
# Build LaTeX table comparing KS distance across scenarios

files <- c(
  small  = "summary_small.csv",
  medium = "summary_medium.csv",
  large  = "summary_large.csv"
)

# Read and tag scenarios
dfs <- lapply(names(files), function(sc) {
  f <- files[[sc]]
  df <- utils::read.csv(f, stringsAsFactors = FALSE)
  df$scenario <- sc
  df
})

df_all <- do.call(rbind, dfs)

# Keep only relevant columns
df_all <- df_all[, c("method", "scenario", "ks")]

# Pivot to wide format (base R)
tab <- reshape(
  df_all,
  idvar = "method",
  timevar = "scenario",
  direction = "wide"
)

# Clean column names
colnames(tab) <- gsub("ks\\.", "", colnames(tab))

# Order methods nicely (optional)
method_order <- c("prior_importance", "importance", "poisson", "power")
tab <- tab[match(method_order, tab$method, nomatch = 0), ]

# Round for presentation
num_cols <- setdiff(colnames(tab), "method")
tab[num_cols] <- lapply(tab[num_cols], function(x) round(x, 3))

# -----------------------------
# LaTeX export
# -----------------------------

latex_lines <- c(
  "\\begin{table}[h]",
  "\\centering",
  "\\begin{tabular}{lccc}",
  "\\hline",
  "Method & Small & Medium & Large \\\\",
  "\\hline"
)

for (i in seq_len(nrow(tab))) {
  row <- tab[i, ]
  line <- sprintf(
    "%s & %.3f & %.3f & %.3f \\\\",
    row$method,
    row$small,
    row$medium,
    row$large
  )
  latex_lines <- c(latex_lines, line)
}

latex_lines <- c(
  latex_lines,
  "\\hline",
  "\\end{tabular}",
  "\\caption{Kolmogorov-Smirnov distance between each estimator and the Stan baseline.}",
  "\\label{tab:ks-comparison}",
  "\\end{table}"
)

# Write to file
writeLines(latex_lines, "ks_table.tex")

cat("LaTeX table written to ks_table.tex\n")
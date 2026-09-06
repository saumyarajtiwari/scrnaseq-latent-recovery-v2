# Item 6: null-calibration for Step 6.6 (Neighborhood Collapse), mean k-NN
# overlap < 0.5. Under "no true relationship between reference and method
# embeddings," a cell's k-NN set in the method's embedding is essentially
# a random k-subset of the other n-1 cells -- overlap with the true k-NN
# set is then a purely combinatorial question (intersection of two random
# k-subsets from n-1 items), depending only on (n_cells, k). No actual
# embedding/matrix simulation needed -- direct resampling of index sets.

set.seed(42)
N_DRAWS <- 200
K_NEIGHBORS <- 15

d <- read.csv("data/processed/step6_6_neighborhood_collapse.csv", stringsAsFactors = FALSE)

n_buckets <- sort(unique(d$n_cells[!is.na(d$n_cells)]))
cat("Calibration buckets:", length(n_buckets), "\n")

calibrate_overlap_null <- function(n, k, n_draws, seed_base) {
  set.seed(seed_base + n)
  pool <- seq_len(n - 1)  # candidate neighbors, excluding self
  replicate(n_draws, {
    mean(sapply(seq_len(min(50, n)), function(i) {
      # Sample 50 "cells" per draw (not all n, for speed) -- overlap
      # distribution doesn't depend on which cell, only on n and k
      true_set <- sample(pool, k)
      random_set <- sample(pool, k)
      length(intersect(true_set, random_set)) / k
    }))
  })
}

calib <- list()
for (n in n_buckets) {
  k_use <- min(K_NEIGHBORS, n - 1)
  null_vals <- calibrate_overlap_null(n, k_use, N_DRAWS, 5000)
  calib[[as.character(n)]] <- list(
    n_cells = n, k_used = k_use,
    null_mean = mean(null_vals), null_median = median(null_vals),
    null_95pct = quantile(null_vals, 0.95), null_99pct = quantile(null_vals, 0.99),
    null_max = max(null_vals)
  )
}
calib_df <- do.call(rbind.data.frame, calib)
cat("\nk-NN overlap null calibration:\n"); print(calib_df)

lookup_95pct <- setNames(calib_df$null_95pct, as.character(calib_df$n_cells))
d$overlap_null_95pct <- lookup_95pct[as.character(d$n_cells)]
# Flag direction is LOW overlap = collapse; "surprisingly non-random" means
# BELOW the null's 95th percentile is not meaningful (that's still within
# chance range) -- the relevant comparison is whether real overlap is
# statistically ABOVE chance (structure preserved) or indistinguishable
# from/below it (collapsed). So the calibrated flag is: real overlap is
# NOT statistically distinguishable from chance (<= null 95th pct).
d$overlap_flag_calibrated <- d$mean_overlap_score <= d$overlap_null_95pct

write.csv(d, "data/processed/step6_6_neighborhood_collapse_calibrated.csv", row.names = FALSE)
write.csv(calib_df, "data/processed/step6_6_overlap_null_calibration.csv", row.names = FALSE)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat("Literal flag (overlap<0.5):", sum(d$overlap_flag, na.rm=TRUE), "\n")
cat("Calibrated flag (overlap indistinguishable from/below chance):", sum(d$overlap_flag_calibrated, na.rm=TRUE), "\n")
cat("\nReal overlap median:", median(d$mean_overlap_score, na.rm=TRUE), "\n")
cat("Null 95th-pct range across buckets:", range(calib_df$null_95pct), "\n")

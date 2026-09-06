# Item 6: null-calibration for Step 6.5a (Over-Smoothing, real-data DPT
# branch), Spearman(DPT, binary reference) < 0.6. Only 18 rows total --
# calibrated exactly per-row (matching each row's actual n_progenitor/
# n_mature split, not just total n, since correlation with a binary
# variable's null variance depends on class balance too) rather than
# bucketing, which would be unnecessary approximation at this scale.

set.seed(42)
N_DRAWS <- 2000  # cheap at this row count -- use more draws for precision

d <- read.csv("data/processed/step6_5a_pancreas_dpt.csv", stringsAsFactors = FALSE)
cat("Total rows:", nrow(d), "\n")

calibrate_null <- function(n_prog, n_mat, n_draws, seed) {
  set.seed(seed)
  ref <- c(rep(0, n_prog), rep(1, n_mat))
  replicate(n_draws, {
    x <- rnorm(n_prog + n_mat)
    suppressWarnings(cor(x, sample(ref), method = "spearman"))
  })
}

d$null_5pct <- NA_real_
d$null_95pct <- NA_real_
for (i in seq_len(nrow(d))) {
  if (is.na(d$n_progenitor[i]) || is.na(d$n_mature[i])) next
  null_vals <- calibrate_null(d$n_progenitor[i], d$n_mature[i], N_DRAWS, 6000 + i)
  d$null_5pct[i] <- quantile(null_vals, 0.05)
  d$null_95pct[i] <- quantile(null_vals, 0.95)
}

# Flag direction: rho < 0.6 = over-smoothing (LOW correlation = bad).
# Calibrated: is real rho below what pure chance would typically produce?
# (i.e., genuinely indistinguishable from/worse than random, not just
# below the literal 0.6 bar)
d$over_smoothing_flag_calibrated <- d$spearman_dpt_vs_reference <= d$null_95pct

write.csv(d, "data/processed/step6_5a_pancreas_dpt_calibrated.csv", row.names = FALSE)

cat("\n=== CALIBRATION COMPLETE ===\n")
print(d[, c("source","method","n_progenitor","n_mature","spearman_dpt_vs_reference",
            "over_smoothing_flag","null_95pct","over_smoothing_flag_calibrated")])
cat("\nLiteral flag (rho<0.6):", sum(d$over_smoothing_flag, na.rm=TRUE), "/", sum(!is.na(d$over_smoothing_flag)), "\n")
cat("Calibrated flag (indistinguishable from/below chance):", sum(d$over_smoothing_flag_calibrated, na.rm=TRUE), "/", sum(!is.na(d$over_smoothing_flag_calibrated)), "\n")

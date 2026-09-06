# Item 6: null-calibration for Step 6.5c (Over-Smoothing, simulated
# variance-ratio branch), ratio < 0.20. This branch found 0/31,640
# flagged under the literal spec -- the informative calibration question
# here is inverted from the other five tests: not "how many flags are
# noise" (there are none), but "is the 0.20 threshold so lenient that
# this test could never detect genuine over-smoothing" -- i.e., does it
# have any statistical power at all. Tested via label-permutation null:
# keep real gene-expression matrix and real embedding fixed, permute
# true_group, recompute both pre_frac/post_frac and their ratio -- this
# is what the ratio would look like if there were NO true biological
# signal. Stratified subset (not full 31,640-row grid, too expensive):
# ~80 files across 4 methods x 3 simulators x 2 target_bands, matching
# Item 3's cross-simulator triangulation standard rather than a single-
# simulator sample.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

set.seed(42)
N_DRAWS <- 200
METHODS <- c("pca_raw", "pca_libnorm", "pca_log", "pca_shiftedlog")
TARGET_SUM <- 10000

signal_fraction <- function(mat, groups, cells_are_rows) {
  if (!cells_are_rows) mat <- t(mat)
  grand_mean <- colMeans(mat)
  total_ss <- sum(sweep(mat, 2, grand_mean, "-")^2)
  if (total_ss == 0) return(NA_real_)
  between_ss <- 0
  for (g in unique(groups)) {
    idx <- groups == g
    n_g <- sum(idx)
    group_mean <- colMeans(mat[idx, , drop = FALSE])
    between_ss <- between_ss + n_g * sum((group_mean - grand_mean)^2)
  }
  between_ss / total_ss
}

derive_embedding_path <- function(raw_path, method) {
  bn <- basename(raw_path)
  sim <- basename(dirname(raw_path))
  suffix_map <- c(pca_raw = "rawpca", pca_libnorm = "libnormpca",
                   pca_log = "logpca", pca_shiftedlog = "shiftedlogpca")
  new_bn <- sub("\\.rds$", paste0("_", suffix_map[[method]], ".rds"), bn)
  file.path("data/processed", method, "simulated", sim, new_bn)
}

d5c <- read.csv("data/processed/step6_5c_simulated_variance_ratio.csv", stringsAsFactors = FALSE)
d5c$sim <- basename(dirname(d5c$file_path))

# Stratified sample: unique raw files, spread across (sim, target_band)
files_meta <- unique(d5c[, c("file_path","sim","target_band")])
strata <- split(files_meta, list(files_meta$sim, files_meta$target_band))
n_per_stratum <- 5  # 3 sims x 2 bands x 5 = 30 files, x4 methods computed each = 120 calibration cases
sample_files <- do.call(rbind, lapply(strata, function(s) {
  if (nrow(s) == 0) return(NULL)
  s[sample(nrow(s), min(n_per_stratum, nrow(s))), ]
}))
cat("Sampled files:", nrow(sample_files), "across strata:\n")
print(table(sample_files$sim, sample_files$target_band))

calibrate_one <- function(raw_path, method, n_draws) {
  src <- readRDS(raw_path)
  counts_mat <- counts(src)
  tg <- as.character(src$true_group)

  lib_sizes <- Matrix::colSums(counts_mat)
  safe_lib <- lib_sizes; safe_lib[safe_lib == 0] <- 1
  scale_factor <- TARGET_SUM / safe_lib
  normalized <- counts_mat %*% Diagonal(x = scale_factor)
  delta_opt <- median(normalized@x)
  S <- normalized
  S@x <- log(S@x + delta_opt) - log(delta_opt)

  pre_pca_transforms <- list(pca_raw = counts_mat, pca_libnorm = normalized,
                              pca_log = log1p(normalized), pca_shiftedlog = S)
  tmat <- as.matrix(pre_pca_transforms[[method]])

  emb_path <- derive_embedding_path(raw_path, method)
  if (!file.exists(emb_path)) return(NULL)
  emb_obj <- readRDS(emb_path)
  embedding <- emb_obj$embedding

  set.seed(7000 + nchar(raw_path) + which(METHODS == method) * 13)
  null_ratios <- replicate(n_draws, {
    perm_tg <- sample(tg)
    pre_frac <- signal_fraction(tmat, perm_tg, cells_are_rows = FALSE)
    post_frac <- signal_fraction(embedding, perm_tg, cells_are_rows = TRUE)
    if (!is.na(pre_frac) && pre_frac > 0) post_frac / pre_frac else NA_real_
  })
  null_ratios[!is.na(null_ratios)]
}

calib_results <- list()
for (i in seq_len(nrow(sample_files))) {
  fp <- sample_files$file_path[i]
  for (method in METHODS) {
    null_vals <- tryCatch(calibrate_one(fp, method, N_DRAWS), error = function(e) NULL)
    if (is.null(null_vals) || length(null_vals) < 10) next
    calib_results[[paste(fp, method)]] <- data.frame(
      file_path = fp, sim = sample_files$sim[i], target_band = sample_files$target_band[i],
      method = method, n_null_draws = length(null_vals),
      null_min = min(null_vals), null_5pct = quantile(null_vals, 0.05),
      null_median = median(null_vals), null_mean = mean(null_vals)
    )
  }
  if (i %% 5 == 0) cat("Progress:", i, "/", nrow(sample_files), "\n")
}
calib_df <- do.call(rbind, calib_results)

write.csv(calib_df, "data/processed/step6_5c_variance_ratio_null_calibration.csv", row.names = FALSE)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat("Calibration cases:", nrow(calib_df), "\n")
cat("\nNull ratio distribution (should be near 1 if the test has power; near/below 0.20 would mean low power):\n")
print(summary(calib_df$null_median))
cat("\nBy method:\n")
print(aggregate(cbind(null_5pct, null_median) ~ method, data = calib_df, FUN = mean))
cat("\nFraction of calibration cases where null 5th pct already falls below 0.20 (would indicate the test could flag pure noise, or conversely that real ratios near 0.20 aren't distinguishable from permuted labels):\n")
cat(mean(calib_df$null_5pct < 0.20, na.rm = TRUE), "\n")

# Item 6: null-calibration for Step 6.7 (Subspace Rotation Slippage), 30 deg
# threshold. Exact closed-form null (no simulation): for a uniformly random
# unit vector in R^p, cos^2(angle to a fixed r-dim subspace) ~ Beta(r/2,
# (p-r)/2) exactly. p (ambient gene dimension) varies per row for
# SCTransform v2/GLM-PCA (read from each file's loadings; found to vary
# substantially below the nominal 3000/2000 cap, likely additional
# zero-count-gene removal at high-sparsity conditions) and is deterministic
# by source for the 4 linear methods (no file reads needed there).

suppressPackageStartupMessages(library(parallel))
ALPHA <- 0.05

d <- read.csv("data/processed/step6_7_subspace_rotation_slippage.csv", stringsAsFactors = FALSE)

# ---- p for the 4 linear methods: deterministic by source ----
p_linear <- c(scdesign3 = 2000, symsim = 2000, splatter = 10000)
d$p_ambient <- NA_integer_
linear_methods <- c("pca_raw","pca_libnorm","pca_log","pca_shiftedlog")
is_linear_sim <- d$method %in% linear_methods & d$source %in% names(p_linear)
d$p_ambient[is_linear_sim] <- p_linear[d$source[is_linear_sim]]

# ---- p for real-data linear-method rows: one representative file read per source ----
real_linear_rows <- d[d$method %in% linear_methods & !(d$source %in% names(p_linear)), ]
if (nrow(real_linear_rows) > 0) {
  real_sources <- unique(real_linear_rows$source)
  for (src in real_sources) {
    fp <- real_linear_rows$file_path[real_linear_rows$source == src][1]
    f <- readRDS(fp)
    p_val <- nrow(f$loadings)
    d$p_ambient[d$source == src & d$method %in% linear_methods] <- p_val
  }
}

# ---- p for SCTransform v2/GLM-PCA: varies per file, read directly, parallelized ----
nl_idx <- which(d$method %in% c("pca_sctransform_v2","pca_glmpca"))
cat("Reading p (ambient gene dim) for", length(nl_idx), "SCTransform/GLM-PCA rows (parallelized)...\n")
t0 <- Sys.time()
p_vals <- mclapply(d$file_path[nl_idx], function(fp) {
  tryCatch(nrow(readRDS(fp)$loadings), error = function(e) NA_integer_)
}, mc.cores = 8)
cat("Elapsed:", format(Sys.time() - t0), "\n")
d$p_ambient[nl_idx] <- unlist(p_vals)

cat("\np_ambient summary:\n"); print(summary(d$p_ambient))
cat("Any NA p_ambient:", sum(is.na(d$p_ambient)), "\n")

# ---- Exact closed-form null threshold per row ----
# theta_null_alpha = angle such that only `alpha` fraction of purely random
# vectors would show an angle this small or smaller (i.e. "surprisingly
# well-aligned by pure chance" cutoff, at the 5% level).
null_angle_threshold <- function(r, p, alpha) {
  if (is.na(r) || is.na(p) || p <= r) return(NA_real_)
  cos2_thresh <- qbeta(1 - alpha, r / 2, (p - r) / 2)
  acos(sqrt(cos2_thresh)) * 180 / pi
}

d$angle_null_5pct_threshold <- mapply(null_angle_threshold, d$rank_true, d$p_ambient, MoreArgs = list(alpha = ALPHA))

# Calibrated flag: observed angle is AT OR BEYOND what pure chance readily
# produces (i.e. NOT surprising given the null -- statistically
# indistinguishable from, or worse than, complete randomness).
d$flag_pc1_calibrated <- d$angle_pc1 >= d$angle_null_5pct_threshold

write.csv(d, "data/processed/step6_7_subspace_rotation_slippage_calibrated.csv", row.names = FALSE)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat("Literal flag (angle_pc1 > 30):", sum(d$flag_pc1, na.rm=TRUE), "\n")
cat("Calibrated flag (angle >= null 5th pct, i.e. indistinguishable from chance):",
    sum(d$flag_pc1_calibrated, na.rm=TRUE), "\n")
cat("\nNull 5th-pct threshold summary (degrees):\n"); print(summary(d$angle_null_5pct_threshold))
cat("\nReal angle_pc1 summary (degrees):\n"); print(summary(d$angle_pc1))

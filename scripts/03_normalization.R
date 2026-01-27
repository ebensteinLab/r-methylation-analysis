#!/usr/bin/env Rscript

# ================================================================
# Script 03: Post-processing / normalization of SeSAMe output
# ================================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

library(matrixStats)

message("Loading SeSAMe-processed matrices...")

mval_matrix <- readRDS("results/processed/mval_matrix_sesame.rds")

# ------------------------------------------------
# OPTIONAL: probe-wise centering / scaling
# (useful for PCA / ML, NOT part of preprocessing)
# ------------------------------------------------

message("Applying optional scaling to M-values...")

# Compute SD per probe
probe_sd <- rowSds(mval_matrix, na.rm = TRUE)

# Count non-NA values per probe
probe_n <- rowCounts(mval_matrix, value = NA, invert = TRUE)

# Keep probes with:
#  - positive variance
#  - at least 2 non-NA values
keep <- probe_sd > 0 & probe_n >= 2

# Sanity check
stopifnot(!any(is.na(keep)))

# Apply scaling
# Subset first
X <- mval_matrix[keep, ]

# Compute row means & SDs
row_means <- rowMeans2(X, na.rm = TRUE)
row_sds   <- rowSds(X, na.rm = TRUE)

# Scale in-place
X <- (X - row_means) / row_sds

mval_scaled <- X
rm(X, row_means, row_sds, mval_matrix, probe_n, probe_sd)
gc()

removed <- sum(!keep)
pct <- round(100 * removed / length(keep), 2)

message("Removed ", removed, " low-information probes (", pct, "%)")

# Check for numerical issues
if (anyNA(mval_scaled)) {
  warning("NAs introduced during scaling; check low-variance probes.")
}

# ------------------------------------------------
# Save outputs
# ------------------------------------------------
saveRDS(mval_scaled, "results/processed/mval_matrix_sesame_scaled.rds")

rm(mval_scaled, keep)
gc()

message("Script 03 completed successfully.")

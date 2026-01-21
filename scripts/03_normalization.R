#!/usr/bin/env Rscript

# ================================================================
# Script 03: Post-processing / normalization of SeSAMe output
# ================================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

message("Loading SeSAMe-processed matrices...")

mval_matrix <- readRDS("results/processed/mval_matrix_sesame.rds")

# ------------------------------------------------
# OPTIONAL: probe-wise centering / scaling
# (useful for PCA / ML, NOT part of preprocessing)
# ------------------------------------------------

message("Applying optional scaling to M-values...")

# Compute SD per probe
probe_sd <- apply(mval_matrix, 1, sd, na.rm = TRUE)

# Count non-NA values per probe
probe_n <- rowSums(!is.na(mval_matrix))

# Keep probes with:
#  - positive variance
#  - at least 2 non-NA values
keep <- probe_sd > 0 & probe_n >= 2

# Sanity check
stopifnot(!any(is.na(keep)))

# Apply scaling
mval_scaled <- t(scale(t(mval_matrix[keep, ])))
rm(mval_matrix)
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

message("Script 03 completed successfully.")

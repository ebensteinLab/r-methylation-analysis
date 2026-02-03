#!/usr/bin/env Rscript

# ============================================================
# 09_cf_fractions_batch_correction.R
# Batch correction of cell fractions using per-cell-type
# linear models
# ============================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(data.table)
})

IN_FILE  <- "results/deconvolution/cell_fractions_idol_nnls.rds"
OUT_DIR  <- "results/deconvolution"
OUT_FILE <- file.path(OUT_DIR, "cell_fractions_idol_nnls_batch_corrected.rds")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading deconvolution results...")
df <- readRDS(IN_FILE)

stopifnot("Patient" %in% colnames(df))
stopifnot("Sentrix_Id" %in% colnames(df))

batch <- factor(df$Sentrix_Id)

# Identify fraction columns
fraction_cols <- setdiff(
  colnames(df),
  c("Patient", "Sentrix_Id", "Disease")
)

message("Cell types: ", paste(fraction_cols, collapse = ", "))
message("Batches: ", nlevels(batch))

corrected <- df

# ------------------------------------------------------------
# Batch correction via linear models (per cell type)
# ------------------------------------------------------------
for (ct in fraction_cols) {
  
  y <- df[[ct]]
  
  # Skip if constant
  if (sd(y, na.rm = TRUE) == 0) {
    message("Skipping ", ct, " (no variance)")
    next
  }
  
  # Linear model
  fit <- lm(y ~ batch)
  
  # Residuals + global mean
  corrected[[ct]] <- residuals(fit) + mean(y, na.rm = TRUE)
}

# ------------------------------------------------------------
# Post-processing: enforce constraints
# ------------------------------------------------------------

# No negatives
corrected[fraction_cols] <- lapply(
  corrected[fraction_cols],
  function(x) pmax(x, 0)
)

# Renormalize to sum = 1
row_sums <- rowSums(corrected[fraction_cols])
row_sums[row_sums == 0] <- 1

corrected[fraction_cols] <- corrected[fraction_cols] / row_sums

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
saveRDS(corrected, OUT_FILE)

csv_out <- sub("\\.rds$", ".csv", OUT_FILE)
fwrite(corrected, csv_out)

message("Saved batch-corrected fractions:")
message(" - ", OUT_FILE)
message(" - ", csv_out)

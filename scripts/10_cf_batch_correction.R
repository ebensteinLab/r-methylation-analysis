#!/usr/bin/env Rscript

# ================================================================
# Script 04: Batch correction using ComBat (SeSAMe pipeline)
# ================================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(sva)
  library(sesame)
  library(GenomicRanges)
})

message("Loading SeSAMe-processed data...")

normalize_sesame_ids <- function(mat) {
  base <- sub("_.*$", "", rownames(mat))
  rownames(mat) <- base
  mat[!duplicated(base), ]
}

cf_mval <- readRDS("results/processed/cf_mval_matrix_raw_sesame.rds")
cf_mval <- normalize_sesame_ids(cf_mval)
cf_targets <- readRDS("results/processed/cf_targets_merged.rds")

message("Number of rows in targets: ", nrow(cf_targets))

# ------------------------------------------------
# Ensure sample alignment
# ------------------------------------------------

cf_mval_samples    <- colnames(cf_mval)
cf_target_samples  <- cf_targets$Patient

missing_in_cf_targets <- setdiff(cf_mval_samples, cf_target_samples)
cf_missing_in_mval    <- setdiff(cf_target_samples, cf_mval_samples)
rm(cf_mval_samples, cf_target_samples)

if (length(missing_in_cf_targets) > 0 || length(cf_missing_in_mval) > 0) {
  
  message("---- Sample name mismatch detected ----")
  
  if (length(missing_in_cf_targets) > 0) {
    message("Samples in mval but NOT in targets (", length(missing_in_cf_targets), "):")
    print(missing_in_cf_targets)
  }
  
  if (length(cf_missing_in_mval) > 0) {
    message("Samples in targets but NOT in mval (", length(cf_missing_in_mval), "):")
    print(cf_missing_in_mval)
  }
  
  stop("Sample names in mval matrix do not match targets table.")
}

cf_targets <- cf_targets[match(colnames(cf_mval), cf_targets$Patient), ]

# ------------------------------------------------
# Define batch variable: cf_targets$Sentrix_Id
# ------------------------------------------------

batch <- as.factor(cf_targets$Sentrix_Id)

if (anyNA(batch)) {
  stop("Batch variable contains NA values.")
}

message("Batch levels: ", paste(levels(batch), collapse = ", "))

keep <- rowMeans(!is.na(cf_mval)) > 0.7
cf_mval <- cf_mval[keep, ]
mod <- model.matrix(~ 1, data = cf_targets)  # no biological covariates

rm(cf_targets, keep)
gc()

# ------------------------------------------------
# Run ComBat
# ------------------------------------------------
if (nlevels(batch) < 2) {
  warning("Only one batch detected; skipping ComBat.")
  mval_corrected <- cf_mval
  combat_ran <- FALSE
} else {
  gr <- sesameData_getManifestGRanges("EPICv2")
  chr <- as.character(seqnames(gr))
  names(chr) <- names(gr)
  
  # Keep only probes in matrix
  chr <- chr[rownames(cf_mval)]
  
  mval_corrected <- matrix(NA, nrow = nrow(cf_mval), ncol = ncol(cf_mval),
                           dimnames = dimnames(cf_mval))
  
  for (c in unique(chr)) {
    message("Running ComBat on ", c)
    
    idx <- which(chr == c)
    X <- cf_mval[idx, ]
    
    # Skip tiny chromosomes
    if (length(idx) < 100) next
    
    Xc <- ComBat(dat = X, batch = batch, mod = mod, par.prior = TRUE)
    mval_corrected[idx, ] <- Xc
    
    rm(X, Xc); gc()
  }
  combat_ran <- TRUE
  rm(gr, chr)
}

rm(cf_mval, mod, batch)
gc()

# ------------------------------------------------
# Save output
# ------------------------------------------------

message("Saving outputs")

mask_matrix <- readRDS("results/processed/cf_mask_matrix_sesame.rds")
mask <- rowMeans(mask_matrix) > 0.7
rm(mask_matrix)
gc()

saveRDS(mval_corrected, "results/processed/cf_mval_matrix_sesame_batch_corrected.rds")

# Convert corrected M-values back to betas
beta_corrected <- MValueToBetaValue(mval_corrected)

rm(mval_corrected, mask)
gc()

saveRDS(beta_corrected, "results/processed/cf_beta_matrix_sesame_batch_corrected.rds")

rm(beta_corrected)
gc()

if (combat_ran) {
  message("Batch correction completed successfully.")
} else {
  message("Batch correction skipped (single batch).")
}

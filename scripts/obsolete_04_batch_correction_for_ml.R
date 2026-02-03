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

mval <- readRDS("results/processed/mval_matrix_raw_sesame.rds")
targets <- readRDS("results/processed/targets_merged.rds")

message("Number of rows in targets: ", nrow(targets))

# ------------------------------------------------
# Ensure sample alignment
# ------------------------------------------------

mval_samples    <- colnames(mval)
target_samples  <- targets$Patient

missing_in_targets <- setdiff(mval_samples, target_samples)
missing_in_mval    <- setdiff(target_samples, mval_samples)

if (length(missing_in_targets) > 0 || length(missing_in_mval) > 0) {
  
  message("---- Sample name mismatch detected ----")
  
  if (length(missing_in_targets) > 0) {
    message("Samples in mval but NOT in targets (", length(missing_in_targets), "):")
    print(missing_in_targets)
  }
  
  if (length(missing_in_mval) > 0) {
    message("Samples in targets but NOT in mval (", length(missing_in_mval), "):")
    print(missing_in_mval)
  }
  
  stop("Sample names in mval matrix do not match targets table.")
}

targets <- targets[match(colnames(mval), targets$Patient), ]

# ------------------------------------------------
# Define batch variable: targets$Sentrix_Id
# ------------------------------------------------

batch <- as.factor(targets$Sentrix_Id)

if (anyNA(batch)) {
  stop("Batch variable contains NA values.")
}

message("Batch levels: ", paste(levels(batch), collapse = ", "))

keep <- rowSums(!is.na(mval)) == length(batch)
mval <- mval[keep, ]
mod <- model.matrix(~ 1, data = targets)  # no biological covariates

rm(targets, keep)
gc()

# ------------------------------------------------
# Run ComBat
# ------------------------------------------------
if (nlevels(batch) < 2) {
  warning("Only one batch detected; skipping ComBat.")
  mval_corrected <- mval
  combat_ran <- FALSE
} else {
  gr <- sesameData_getManifestGRanges("EPICv2")
  chr <- as.character(seqnames(gr))
  names(chr) <- names(gr)
  
  # Keep only probes in matrix
  chr <- chr[rownames(mval)]
  
  mval_corrected <- matrix(NA, nrow = nrow(mval), ncol = ncol(mval),
                           dimnames = dimnames(mval))
  
  for (c in unique(chr)) {
    message("Running ComBat on ", c)
    
    idx <- which(chr == c)
    X <- mval[idx, ]
    
    # Skip tiny chromosomes
    if (length(idx) < 100) next
    
    Xc <- ComBat(dat = X, batch = batch, mod = mod, par.prior = TRUE)
    mval_corrected[idx, ] <- Xc
    
    rm(X, Xc); gc()
  }
  combat_ran <- TRUE
  rm(gr, chr, idx)
}

rm(mval, mod, batch)
gc()

# ------------------------------------------------
# Save output
# ------------------------------------------------

message("Saving outputs")

mask_full <- readRDS("results/processed/mask_matrix_sesame.rds")
mask <- mask_full[rownames(mval_corrected), , drop = FALSE]
mval_corrected[mask == 0] <- NA

saveRDS(mval_corrected, "results/processed/mval_matrix_sesame_batch_corrected.rds")

# Convert corrected M-values back to betas
beta_corrected <- MValueToBetaValue(mval_corrected)

rm(mval_corrected, mask, mask_full)
gc()

saveRDS(beta_corrected, "results/processed/beta_matrix_sesame_batch_corrected.rds")

rm(beta_corrected)
gc()

if (combat_ran) {
  message("Batch correction completed successfully.")
} else {
  message("Batch correction skipped (single batch).")
}

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
})

message("Loading SeSAMe-processed data...")

mval <- readRDS("results/processed/mval_matrix_raw_sesame.rds")
targets <- readRDS("results/processed/targets_merged.rds")

message("Number of rows in targets: ", nrow(targets))

# ------------------------------------------------
# Ensure sample alignment
# ------------------------------------------------

if (!all(colnames(mval) %in% targets$sample_name)) {
  stop("Sample names in mval matrix do not match targets table.")
}

targets <- targets[match(colnames(mval), targets$sample_name), ]

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
  message("Running ComBat with ", nlevels(batch), " batches.")
  mval_corrected <- ComBat(
    dat = mval,
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  combat_ran <- TRUE
}

rm(mval, mod, batch)
gc()

# ------------------------------------------------
# Save output
# ------------------------------------------------

mask <- readRDS("results/processed/final_probe_mask.rds")
mval_corrected[!mask, ] <- NA

saveRDS(mval_corrected, "results/processed/mval_matrix_sesame_batch_corrected.rds")

# Convert corrected M-values back to betas
beta_corrected <- MValueToBetaValue(mval_corrected)

rm(mval_corrected)
gc()

saveRDS(beta_corrected, "results/processed/beta_matrix_sesame_batch_corrected.rds")

rm(beta_corrected)
gc()

if (combat_ran) {
  message("Batch correction completed successfully.")
} else {
  message("Batch correction skipped (single batch).")
}

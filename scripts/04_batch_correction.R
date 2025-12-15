#!/usr/bin/env Rscript

# ================================================================
# Script 04: Batch correction using ComBat (SeSAMe pipeline)
# ================================================================

suppressPackageStartupMessages({
  library(sva)
})

message("Loading SeSAMe-processed data...")

mval <- readRDS("results/processed/mval_matrix_sesame_scaled.rds")
targets <- readRDS("results/processed/targets_with_sesame.rds")

# ------------------------------------------------
# Ensure sample alignment
# ------------------------------------------------

if (!all(colnames(mval) %in% targets$sample_name)) {
  stop("Sample names in mval matrix do not match targets table.")
}

targets <- targets[match(colnames(mval), targets$sample_name), ]

# ------------------------------------------------
# Define batch variable
# ------------------------------------------------
# Example batch variables:
#   targets$Batch
#   targets$Sentrix_Id
#   targets$RunDate
#   targets$Study

batch <- as.factor(targets$Batch)

if (anyNA(batch)) {
  stop("Batch variable contains NA values.")
}

message("Batch levels: ", paste(levels(batch), collapse = ", "))

# ------------------------------------------------
# OPTIONAL: model biological covariates
# ------------------------------------------------
# Example:
# mod <- model.matrix(~ Disease + Age + Sex, data = targets)

mod <- model.matrix(~ 1, data = targets)  # no biological covariates

# ------------------------------------------------
# Run ComBat
# ------------------------------------------------

batch <- as.factor(targets$Batch)

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

# ------------------------------------------------
# Save output
# ------------------------------------------------

saveRDS(
  mval_corrected,
  "results/processed/mval_matrix_sesame_batch_corrected.rds"
)

# Convert corrected M-values back to betas
beta_corrected <- MValueToBetaValue(mval_corrected)

saveRDS(
  beta_corrected,
  "results/processed/beta_matrix_sesame_batch_corrected.rds"
)

if (combat_ran) {
  message("Batch correction completed successfully.")
} else {
  message("Batch correction skipped (single batch).")
}

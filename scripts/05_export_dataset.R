#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(data.table)
})

message("Loading batch-corrected methylation matrices...")

# ------------------------------------------------
# Load data
# ------------------------------------------------
mval_corrected <- readRDS(
  "results/processed/mval_matrix_sesame_batch_corrected.rds"
)

beta_corrected <- readRDS("results/processed/beta_matrix_sesame_batch_corrected.rds")

targets <- readRDS("results/processed/targets_merged.rds")

# ------------------------------------------------
# Sanity checks
# ------------------------------------------------
stopifnot(
  identical(colnames(mval_corrected), colnames(beta_corrected)),
  all(colnames(mval_corrected) %in% targets$sample_name)
)

# Reorder metadata to match matrix columns
targets <- targets[
  match(colnames(mval_corrected), targets$sample_name),
]

# Final check
stopifnot(
  identical(colnames(mval_corrected), targets$sample_name)
)

# ------------------------------------------------
# Export CSVs
# ------------------------------------------------
message("Writing output files...")

fwrite(
  as.data.table(mval_corrected, keep.rownames = "Probe_ID"),
  "results/processed/mval_matrix_batch_corrected.csv"
)

fwrite(
  as.data.table(beta_corrected, keep.rownames = "Probe_ID"),
  "results/processed/beta_matrix_batch_corrected.csv"
)

fwrite(targets, "results/processed/sample_metadata.csv")

rm(beta_corrected, mval_corrected, targets)
gc()

message("Export completed successfully.")

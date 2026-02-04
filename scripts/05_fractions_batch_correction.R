#!/usr/bin/env Rscript

# ================================================================
# Script 05: Batch correction of IDOL deconvolution fractions
# ================================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(sva)
  library(compositions)
})

message("Loading deconvolution results...")

epsilon <- 1e-6

# ------------------------------------------------
# Load inputs
# ------------------------------------------------
fractions_df <- readRDS("results/deconvolution/blood_cell_fractions_idol_genomic.rds")

targets <- readRDS("results/processed/targets_merged.rds")

# ------------------------------------------------
# Sanity checks
# ------------------------------------------------
stopifnot(
  "Patient" %in% colnames(fractions_df),
  "Sentrix_Id" %in% colnames(targets)
)

# Align targets to fractions
targets <- targets[match(fractions_df$Patient, targets$Patient),]

stopifnot(identical(fractions_df$Patient, targets$Patient))

# ------------------------------------------------
# Extract fraction matrix
# ------------------------------------------------
fraction_cols <- setdiff(
  colnames(fractions_df),
  c("Patient", "Disease")
)

fractions <- as.matrix(fractions_df[, fraction_cols])
rownames(fractions) <- fractions_df$Patient

stopifnot(is.numeric(fractions))

# ------------------------------------------------
# Define batch variable
# ------------------------------------------------
batch <- as.factor(targets$Sentrix_Id)

if (anyNA(batch)) {
  stop("Batch variable (Sentrix_Id) contains NA values.")
}

message("Batch levels: ", paste(levels(batch), collapse = ", "))

if (nlevels(batch) < 2) {
  warning("Only one batch detected; skipping ComBat.")
  fractions_corrected <- fractions
  combat_ran <- FALSE
} else {
  # No biological covariates here
  mod <- model.matrix(~ 1, data = targets)
  
  fractions_clr <- clr(fractions + epsilon)
  
  message("Running ComBat on fraction matrix...")
  fractions_clr_corrected <- ComBat(
    dat = t(fractions_clr),
    batch = batch,
    mod = mod,
    par.prior = TRUE
  )
  
  fractions_clr_corrected <- t(fractions_clr_corrected)
  fractions_corrected <- exp(fractions_clr_corrected)
  fractions_corrected <- fractions_corrected / rowSums(fractions_corrected)
  
  combat_ran <- TRUE
}

# ------------------------------------------------
# Reassemble output
# ------------------------------------------------
fractions_corrected_df <- cbind(
  fractions_df[, c("Patient"), drop = FALSE],
  as.data.frame(fractions_corrected)
)

if ("Disease" %in% colnames(fractions_df)) {
  fractions_corrected_df$Disease <- fractions_df$Disease
}

# ------------------------------------------------
# Save outputs
# ------------------------------------------------
out_rds <- "results/deconvolution/blood_cell_fractions_idol_genomic_batch_corrected.rds"
out_csv <- "results/deconvolution/blood_cell_fractions_idol_genomic_batch_corrected.csv"

saveRDS(fractions_corrected_df, out_rds)

# ------------------------------------------------
# Round numeric columns for CSV output (5 digits)
# ------------------------------------------------
fractions_corrected_csv <- fractions_corrected_df

num_cols <- sapply(fractions_corrected_csv, is.numeric)
fractions_corrected_csv[, num_cols] <-
  round(fractions_corrected_csv[, num_cols], 5)
write.csv(fractions_corrected_csv, out_csv, row.names = FALSE)

message("Saved batch-corrected fractions to:")
message(out_rds)
message(out_csv)

if (combat_ran) {
  message("Batch correction completed successfully.")
} else {
  message("Batch correction skipped (single batch).")
}

rm(fractions, fractions_corrected, fractions_df, targets)
gc()

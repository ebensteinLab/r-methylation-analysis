#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(EpiDISH)
})

message("Loading inputs...")

# ------------------------------------------------
# Load bulk beta values and metadata
# ------------------------------------------------
beta <- readRDS("results/processed/beta_matrix_sesame_batch_corrected.rds")
targets <- readRDS("results/processed/targets_merged.rds")

stopifnot(ncol(beta) == nrow(targets))

# ------------------------------------------------
# Load IDOL outputs
# ------------------------------------------------
idol_probes <- readRDS("results/processed/idol_cpgs.rds")
epicv2_to_epic <- readRDS("results/processed/epicv2_to_epic_map.rds")
ref_beta_centroids <- readRDS("results/processed/ref_beta_epic_centroids.rds")

message("IDOL CpGs: ", length(idol_probes))
message("Reference cell types: ", paste(colnames(ref_beta_centroids), collapse = ", "))

# ------------------------------------------------
# Harmonize probe IDs (EPICv2 → EPIC v1)
# ------------------------------------------------
old_ids <- rownames(beta)
new_ids <- epicv2_to_epic[old_ids]

keep <- !is.na(new_ids)

beta <- beta[keep, , drop = FALSE]
rownames(beta) <- new_ids[keep]

rm(old_ids, new_ids, keep, epicv2_to_epic)
gc()

# ------------------------------------------------
# Restrict to IDOL CpGs
# ------------------------------------------------
common_idol <- Reduce(
  intersect,
  list(
    idol_probes,
    rownames(beta),
    rownames(ref_beta_centroids)
  )
)

message("CpGs used for deconvolution: ", length(common_idol))
stopifnot(length(common_idol) > 200)

beta_idol <- beta[common_idol, , drop = FALSE]
ref_idol  <- ref_beta_centroids[common_idol, , drop = FALSE]

rm(beta, idol_probes, ref_beta_centroids, common_idol)
gc()

# ------------------------------------------------
# Run EpiDISH (RPC)
# ------------------------------------------------
message("Running EpiDISH (RPC)...")

deconv <- epidish(
  beta.m = beta_idol,
  ref.m  = ref_idol,
  method = "RPC"
)

fractions <- deconv$estF

rm(beta_idol, ref_idol, deconv)
gc()

# ------------------------------------------------
# Assemble final output table
# ------------------------------------------------
fractions_df <- cbind(
  targets[, c("sample_name", "Patient"), drop = FALSE],
  fractions
)

rm(fractions, targets)
gc()

# Optional: derive Disease from Patient (ignore first 5 chars)
fractions_df$Disease <- substr(fractions_df$Patient, 6, nchar(fractions_df$Patient))

# ------------------------------------------------
# Save
# ------------------------------------------------
out_file <- "results/processed/blood_cell_fractions_idol_all_samples.rds"
saveRDS(fractions_df, out_file)

fractions_df_rounded <- fractions_df
num_cols <- sapply(fractions_df_rounded, is.numeric)
fractions_df_rounded[, num_cols] <- round(fractions_df_rounded[, num_cols], 5)

csv_file <- "results/processed/blood_cell_fractions_idol_all_samples.csv"
write.csv(fractions_df_rounded, file = csv_file, row.names = FALSE)

message("Saved deconvolution results to:")
message(out_file)
message(csv_file)

rm(fractions_df, fractions_df_rounded, num_cols, out_file, csv_file)
gc()

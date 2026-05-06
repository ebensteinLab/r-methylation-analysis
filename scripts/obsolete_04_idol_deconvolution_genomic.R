#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(EpiDISH)
})

message("Loading inputs...")

MAX_NA_FRAC <- 0.10   # allow up to 10% missing IDOL probes

# ------------------------------------------------
# Load bulk beta values (should be EPICv2 probe IDs, incl suffix)
# ------------------------------------------------
beta <- readRDS("results/processed/beta_matrix_sesame.rds")
targets <- readRDS("results/processed/targets_with_sesame.rds")

stopifnot(ncol(beta) == nrow(targets))

message("Bulk beta dims: ", nrow(beta), " x ", ncol(beta))
message("Bulk beta rowname example: ", rownames(beta)[1])

# ------------------------------------------------
# Load IDOL outputs (EPICv2 space)
# ------------------------------------------------
idol_probes <- readRDS("results/processed/idol_cpgs_epicv2.rds")
ref_beta_centroids <- readRDS("results/processed/ref_beta_epicv2_centroids.rds")

message("IDOL probes (EPICv2): ", length(idol_probes))
message("Reference cell types: ", paste(colnames(ref_beta_centroids), collapse = ", "))

# ------------------------------------------------
# Restrict to IDOL probes present in BOTH bulk and reference
# ------------------------------------------------
common_idol <- Reduce(
  intersect,
  list(
    idol_probes,
    rownames(beta),
    rownames(ref_beta_centroids)
  )
)

message("Probes used for deconvolution: ", length(common_idol))

# For safety; you can lower temporarily if you’re debugging
stopifnot(length(common_idol) > 200)

beta_idol <- beta[common_idol, , drop = FALSE]

# ------------------------------------------------
# NA diagnostics per sample
# ------------------------------------------------
na_per_sample <- colSums(is.na(beta_idol))
na_frac_per_sample <- na_per_sample / nrow(beta_idol)

summary(na_frac_per_sample)

# Inspect worst samples
print(head(sort(na_frac_per_sample, decreasing = TRUE), 20))

keep_samples <- na_frac_per_sample <= MAX_NA_FRAC

message("Keeping ", sum(keep_samples), " / ", length(keep_samples), " samples")
message("Dropping ", sum(!keep_samples), " samples due to high NA fraction")

dropped_samples <- names(keep_samples)[!keep_samples]
if (length(dropped_samples) > 0) {
  message("Dropped samples:")
  print(dropped_samples)
}

write.csv(data.frame(Patient = dropped_samples), "results/deconvolution/excluded_samples_high_na.csv", row.names = FALSE)

beta_idol <- beta_idol[, keep_samples, drop = FALSE]
targets   <- targets[match(colnames(beta_idol), targets$Patient), ]

stopifnot(identical(colnames(beta_idol), targets$Patient))

ref_idol  <- ref_beta_centroids[common_idol, , drop = FALSE]

rm(beta, idol_probes, ref_beta_centroids, common_idol, dropped_samples)
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
# Assemble output table
# ------------------------------------------------
# Ensure rows of fractions match samples
stopifnot(identical(rownames(fractions), targets$Patient) ||
            identical(rownames(fractions), rownames(targets)) ||
            nrow(fractions) == nrow(targets))

fractions_df <- cbind(
  targets[, c("Patient"), drop = FALSE],
  fractions
)

# Optional: derive Disease from Patient (ignore first 5 chars)
fractions_df$Disease <- substr(fractions_df$Patient, 6, pmin(nchar(fractions_df$Patient), 8))

# ------------------------------------------------
# Save
# ------------------------------------------------
out_file <- "results/deconvolution/blood_cell_fractions_idol_genomic.rds"
saveRDS(fractions_df, out_file)

fractions_df_rounded <- fractions_df
num_cols <- sapply(fractions_df_rounded, is.numeric)
fractions_df_rounded[, num_cols] <- round(fractions_df_rounded[, num_cols], 5)

csv_file <- "results/deconvolution/blood_cell_fractions_idol_genomic.csv"
write.csv(fractions_df_rounded, file = csv_file, row.names = FALSE)

message("Saved deconvolution results to:")
message(out_file)
message(csv_file)

rm(fractions_df, fractions_df_rounded, num_cols)
gc()

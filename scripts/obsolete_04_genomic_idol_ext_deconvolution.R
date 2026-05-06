#!/usr/bin/env Rscript

# ================================================================
# Script 04: Genomic deconvolution of batch-corrected data
#            using published EPIC IDOL-Ext
#
# Input:
#   - results/processed/beta_matrix_sesame_batch_corrected.rds
#   - results/processed/targets_with_sesame.rds
#
# Output:
#   - results/deconvolution/blood_cell_fractions_idol_ext_genomic_batch_corrected_12.rds
#   - results/deconvolution/blood_cell_fractions_idol_ext_genomic_batch_corrected_12.csv
#   - results/deconvolution/excluded_samples_high_na_idol_ext_12.csv
#
# Notes:
#   - Uses published EPIC IDOL-Ext probes and centroid matrix
#   - Harmonizes EPICv2 probe IDs with suffixes to plain cg IDs
#   - Runs EpiDISH RPC deconvolution
# ================================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(ExperimentHub)
  library(BiocFileCache)
  library(EpiDISH)
})

message("Loading inputs...")

MAX_NA_FRAC <- 0.20   # allow up to 20% missing EPIC IDOL-Ext probes

OUT_DECONV_DIR <- "results/deconvolution"
dir.create(OUT_DECONV_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helper: cg00000109_TC21 -> cg00000109
# ------------------------------------------------
normalize_sesame_ids <- function(x) {
  sub("_.*$", "", x)
}

# ------------------------------------------------
# Loading IDOL-Ext reference data (Robust CDN Method)
# ------------------------------------------------
message("Fetching reference data directly from Bioconductor CDN...")

# Explicit URLs for the rda/rds objects
probes_url <- "https://experimenthub.bioconductor.org/fetch/1136"
matrix_url <- "https://experimenthub.bioconductor.org/fetch/1137"

# Create safe temp paths
tmp_p <- tempfile(fileext = ".rds")
tmp_m <- tempfile(fileext = ".rds")

# Increase timeout for potential slow connections
options(timeout = 600)

# Download using libcurl which is more robust for binary files on Linux
download.file(probes_url, tmp_p, method = "libcurl", mode = "wb", quiet = TRUE)
download.file(matrix_url, tmp_m, method = "libcurl", mode = "wb", quiet = TRUE)

# Try loading. If readRDS fails, we try load() as a fallback 
# (in case the Hub sends a wrapped .RData/rda instead of a raw rds)
load_bioc_resource <- function(path) {
  out <- try(readRDS(path), silent = TRUE)
  if (inherits(out, "try-error")) {
    # Fallback for .rda format
    env <- new.env()
    load(path, envir = env)
    out <- as.list(env)[[1]]
  }
  return(out)
}

idol_probes_cg <- load_bioc_resource(tmp_p)
ref_beta_centroids_cg <- load_bioc_resource(tmp_m)

# Check the structure of what was loaded
message("Class: ", class(ref_beta_centroids_cg))
message("Dimensions: ", paste(dim(ref_beta_centroids_cg), collapse = " x "))

# If it's a list, look at the first few elements
if (is.list(ref_beta_centroids_cg) && !is.data.frame(ref_beta_centroids_cg)) {
  print(names(ref_beta_centroids_cg))
}

# Peek at the actual data
print(head(ref_beta_centroids_cg))

# Clean up
unlink(c(tmp_p, tmp_m))

if (is.null(idol_probes_cg) || is.null(ref_beta_centroids_cg)) {
  stop("Failed to parse reference data from CDN.")
}

message("Successfully loaded ", length(idol_probes_cg), " IDOL probes.")
message("Reference cell types: ", paste(colnames(ref_beta_centroids_cg), collapse = ", "))

# ------------------------------------------------
# Load batch-corrected genomic beta values
# ------------------------------------------------
beta <- readRDS("results/processed/beta_matrix_sesame_batch_corrected.rds")
targets <- readRDS("results/processed/targets_with_sesame.rds")

stopifnot(ncol(beta) == nrow(targets))

message("Batch-corrected beta dims: ", nrow(beta), " x ", ncol(beta))
message("Bulk beta rowname example: ", rownames(beta)[1])

# Align targets to beta columns using Patient
if (!"Patient" %in% colnames(targets)) {
  stop("targets must contain a 'Patient' column")
}

targets$Patient <- as.character(targets$Patient)
targets <- targets[match(colnames(beta), targets$Patient), , drop = FALSE]

if (!identical(colnames(beta), targets$Patient)) {
  stop("Column names of beta do not match targets$Patient after alignment")
}

# ------------------------------------------------
# Harmonize beta rownames to cg IDs
# ------------------------------------------------
beta_cg <- normalize_sesame_ids(rownames(beta))

# Map plain cg ID -> first EPICv2 full probe ID present in bulk
cg_to_full <- tapply(rownames(beta), beta_cg, function(v) v[[1]])

# ------------------------------------------------
# Restrict to probes present in:
#   1) published EPIC IDOL-Ext library
#   2) batch-corrected bulk matrix
#   3) published reference centroids
# ------------------------------------------------
common_cg <- Reduce(
  intersect,
  list(
    idol_probes_cg,
    beta_cg,
    rownames(ref_beta_centroids_cg)
  )
)

message("Shared EPIC IDOL-Ext probes (cg IDs): ", length(common_cg))

if (length(common_cg) < 100) {
  stop("Too few shared EPIC IDOL-Ext probes: ", length(common_cg))
}

# Convert shared cg IDs back to EPICv2-style IDs used by bulk matrix
common_epicv2 <- unname(cg_to_full[common_cg])
keep_probe_map <- !is.na(common_epicv2)

common_cg <- common_cg[keep_probe_map]
common_epicv2 <- common_epicv2[keep_probe_map]

message("Shared EPIC IDOL-Ext probes mapped to EPICv2 IDs: ", length(common_epicv2))

if (length(common_epicv2) < 100) {
  stop("Too few mapped EPICv2 probes after cg -> EPICv2 conversion: ", length(common_epicv2))
}

# ------------------------------------------------
# Build beta and reference matrices with matched rownames/order
# ------------------------------------------------
beta_idol <- beta[common_epicv2, , drop = FALSE]
ref_idol <- ref_beta_centroids_cg[common_cg, , drop = FALSE]
rownames(ref_idol) <- common_epicv2

stopifnot(identical(rownames(beta_idol), rownames(ref_idol)))

# ------------------------------------------------
# NA diagnostics per sample
# ------------------------------------------------
na_per_sample <- colSums(is.na(beta_idol))
na_frac_per_sample <- na_per_sample / nrow(beta_idol)

message("NA fraction summary across samples:")
print(summary(na_frac_per_sample))

message("Worst samples by NA fraction:")
print(head(sort(na_frac_per_sample, decreasing = TRUE), 20))

keep_samples <- na_frac_per_sample <= MAX_NA_FRAC

message("Keeping ", sum(keep_samples), " / ", length(keep_samples), " samples")
message("Dropping ", sum(!keep_samples), " samples due to high NA fraction")

dropped_samples <- names(keep_samples)[!keep_samples]
if (length(dropped_samples) > 0) {
  message("Dropped samples:")
  print(dropped_samples)
}

write.csv(
  data.frame(Patient = dropped_samples),
  file.path(OUT_DECONV_DIR, "excluded_samples_high_na_idol_ext_12.csv"),
  row.names = FALSE
)

beta_idol <- beta_idol[, keep_samples, drop = FALSE]
targets <- targets[match(colnames(beta_idol), targets$Patient), , drop = FALSE]

stopifnot(identical(colnames(beta_idol), targets$Patient))

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

# ------------------------------------------------
# Assemble output table
# ------------------------------------------------
stopifnot(
  identical(rownames(fractions), targets$Patient) ||
    identical(rownames(fractions), rownames(targets)) ||
    nrow(fractions) == nrow(targets)
)

fractions_df <- cbind(
  targets[, "Patient", drop = FALSE],
  fractions
)

# Derive Disease if not present
if ("Disease" %in% colnames(targets)) {
  fractions_df$Disease <- targets$Disease
} else {
  fractions_df$Disease <- substr(
    fractions_df$Patient,
    6,
    pmin(nchar(fractions_df$Patient), 8)
  )
}

# ------------------------------------------------
# Save
# ------------------------------------------------
out_file <- file.path(
  OUT_DECONV_DIR,
  "blood_cell_fractions_idol_ext_genomic_batch_corrected_12.rds"
)
saveRDS(fractions_df, out_file)

fractions_df_rounded <- fractions_df
num_cols <- sapply(fractions_df_rounded, is.numeric)
fractions_df_rounded[, num_cols] <- round(fractions_df_rounded[, num_cols], 5)

csv_file <- file.path(
  OUT_DECONV_DIR,
  "blood_cell_fractions_idol_ext_genomic_batch_corrected_12.csv"
)
write.csv(fractions_df_rounded, file = csv_file, row.names = FALSE)

message("Saved deconvolution results to:")
message(out_file)
message(csv_file)

rm(
  beta, beta_cg, beta_idol,
  ref_beta_centroids_cg, ref_idol,
  idol_probes_cg, common_cg, common_epicv2,
  cg_to_full, keep_probe_map,
  deconv, fractions, fractions_df, fractions_df_rounded,
  na_per_sample, na_frac_per_sample, keep_samples,
  dropped_samples, num_cols
)
gc()

message("EPIC IDOL-Ext deconvolution of batch-corrected genomic data completed successfully.")

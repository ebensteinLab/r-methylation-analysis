#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(ExperimentHub)
  library(BiocFileCache)
  library(EpiDISH)
})

message("Loading inputs...")

MAX_NA_FRAC <- 0.40   # allow up to 40% missing probes

OUT_DECONV_DIR <- "results/deconvolution"
dir.create(OUT_DECONV_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helper: cg00000109_TC21 -> cg00000109
# ------------------------------------------------
normalize_sesame_ids <- function(x) {
  sub("_.*$", "", x)
}

# ------------------------------------------------
# Load custom reference panel from CSV
# ------------------------------------------------
message("Loading custom reference panel...")

probe_library <- "assaf"
ref_df <- read.csv(sprintf("gss/%s_reference_matrix.csv", probe_library), check.names = FALSE)

rownames(ref_df) <- ref_df$CpG
ref_df$CpG <- NULL

ref_beta_centroids_cg <- as.matrix(ref_df)
mode(ref_beta_centroids_cg) <- "numeric"

idol_probes_cg <- rownames(ref_beta_centroids_cg)

if (is.null(idol_probes_cg) || length(idol_probes_cg) == 0) {
  stop("Custom reference matrix must have CpG IDs as rownames or a CpG column.")
}

message("Successfully loaded ", length(idol_probes_cg), " custom reference probes.")
message("Reference cell types: ", paste(colnames(ref_beta_centroids_cg), collapse = ", "))

# Check the structure of what was loaded
message("Class: ", class(ref_beta_centroids_cg))
message("Dimensions: ", paste(dim(ref_beta_centroids_cg), collapse = " x "))

# If it's a list, look at the first few elements
if (is.list(ref_beta_centroids_cg) && !is.data.frame(ref_beta_centroids_cg)) {
  print(names(ref_beta_centroids_cg))
}

# Peek at the actual data
print(head(ref_beta_centroids_cg))

if (is.null(idol_probes_cg) || is.null(ref_beta_centroids_cg)) {
  stop("Failed to parse reference data from CDN.")
}

message("Successfully loaded ", length(idol_probes_cg), " IDOL probes.")
message("Reference cell types: ", paste(colnames(ref_beta_centroids_cg), collapse = ", "))

# ------------------------------------------------
# Load batch-corrected genomic beta values
# ------------------------------------------------
beta <- readRDS("results/processed/cf_beta_matrix_sesame_batch_corrected.rds")
targets <- readRDS("results/processed/cf_targets_with_sesame.rds")

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

message("Number of unique cg IDs in cf beta: ", length(unique(beta_cg)))
message("Number of reference cg IDs: ", length(rownames(ref_beta_centroids_cg)))

shared_cg <- intersect(rownames(ref_beta_centroids_cg), unique(beta_cg))
message("Shared cg IDs before mapping back to EPICv2 IDs: ", length(shared_cg))

missing_ref <- setdiff(rownames(ref_beta_centroids_cg), unique(beta_cg))
message("Reference probes missing from cf beta: ", length(missing_ref))
print(head(missing_ref, 20))

# Map plain cg ID -> first EPICv2 full probe ID present in bulk
cg_to_full <- tapply(rownames(beta), beta_cg, function(v) v[[1]])

# ------------------------------------------------
# Find shared probes, preserving reference-panel order
# ------------------------------------------------
common_cg <- rownames(ref_beta_centroids_cg)[rownames(ref_beta_centroids_cg) %in% beta_cg]

message("Shared reference probes (cg IDs): ", length(common_cg))

if (length(common_cg) < 100) {
  stop("Too few shared probes: ", length(common_cg))
}

# Map reference cg IDs to EPICv2 full probe IDs in bulk matrix
common_epicv2 <- unname(cg_to_full[common_cg])

keep_probe_map <- !is.na(common_epicv2)
common_cg <- common_cg[keep_probe_map]
common_epicv2 <- common_epicv2[keep_probe_map]

message("Shared probes mapped to EPICv2 IDs: ", length(common_epicv2))

if (length(common_epicv2) < 100) {
  stop("Too few mapped EPICv2 probes after cg -> EPICv2 conversion: ", length(common_epicv2))
}

# ------------------------------------------------
# Build beta and reference matrices with matched order
# ------------------------------------------------
beta_idol <- beta[common_epicv2, , drop = FALSE]
ref_idol  <- ref_beta_centroids_cg[common_cg, , drop = FALSE]

# Rename reference rows to match the EPICv2 IDs in bulk
rownames(ref_idol) <- common_epicv2

# Final explicit reordering safeguard
ref_idol <- ref_idol[rownames(beta_idol), , drop = FALSE]

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
  file.path(OUT_DECONV_DIR, sprintf("excluded_samples_cf_high_na_%s_12.csv", probe_library)),
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
  sprintf("blood_cell_fractions_%s_cf_batch_corrected_12.rds", probe_library)
)
saveRDS(fractions_df, out_file)

fractions_df_rounded <- fractions_df
num_cols <- sapply(fractions_df_rounded, is.numeric)
fractions_df_rounded[, num_cols] <- round(fractions_df_rounded[, num_cols], 5)

csv_file <- file.path(
  OUT_DECONV_DIR,
  sprintf("blood_cell_fractions_%s_cf_batch_corrected_12.csv", probe_library)
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

message("EPIC deconvolution of batch-corrected cf data completed successfully.")

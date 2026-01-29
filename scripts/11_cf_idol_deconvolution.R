#!/usr/bin/env Rscript

# ============================================================
# 11_cf_idol_deconvolution.R
# Deconvolve (gDNA + cfDNA) bulk EPICv2 SeSAMe matrices using:
# - IDOL probe library (idol_cpgs_shared.rds)
# - Reference centroids on shared CpGs (ref_beta_shared_centroids.rds)
#
# Key fixes:
# - Normalize SeSAMe suffix IDs
# - Collapse duplicates after normalization
# - Restrict to IDOL probes that exist in your bulk matrix
# - Robust NNLS with sum-to-one renorm
# ============================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(data.table)
})

OUT_PROC_DIR <- "results/processed"
OUT_RES_DIR  <- "results/deconvolution"
dir.create(OUT_RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

normalize_sesame_ids_vec <- function(x) sub("_.*$", "", x)

collapse_duplicate_rows_mean <- function(mat) {
  stopifnot(!is.null(rownames(mat)))
  rn <- rownames(mat)
  if (!anyDuplicated(rn)) return(mat)
  counts <- table(rn)
  summed <- rowsum(mat, group = rn, reorder = FALSE)
  means <- summed / as.numeric(counts[rownames(summed)])
  means
}

# Solve non-negative least squares and renormalize to sum=1
solve_fractions_nnls <- function(Y, R) {
  # Y: probes x samples
  # R: probes x celltypes (centroids)
  stopifnot(nrow(Y) == nrow(R))
  
  # Prefer nnls if available
  if (requireNamespace("nnls", quietly = TRUE)) {
    X <- matrix(NA_real_, nrow = ncol(R), ncol = ncol(Y))
    rownames(X) <- colnames(R)
    colnames(X) <- colnames(Y)
    
    Rt <- as.matrix(R)
    for (j in seq_len(ncol(Y))) {
      fit <- nnls::nnls(Rt, Y[, j])
      x <- as.numeric(fit$x)
      x[x < 0] <- 0
      s <- sum(x)
      if (s > 0) x <- x / s
      X[, j] <- x
    }
    return(X)
  }
  
  # Fallback: unconstrained LS -> truncate -> renorm
  warning("Package 'nnls' not found. Falling back to nonnegative LS + renormalization.")
  X <- matrix(NA_real_, nrow = ncol(R), ncol = ncol(Y))
  rownames(X) <- colnames(R)
  colnames(X) <- colnames(Y)
  
  Rt <- as.matrix(R)
  for (j in seq_len(ncol(Y))) {
    x <- tryCatch({
      coef(lm.fit(x = Rt, y = Y[, j]))
    }, error = function(e) rep(NA_real_, ncol(Rt)))
    
    x[is.na(x)] <- 0
    x[x < 0] <- 0
    s <- sum(x)
    if (s > 0) x <- x / s
    X[, j] <- x
  }
  X
}

# ------------------------------------------------------------
# 1) Load inputs
# ------------------------------------------------------------
idol_file <- file.path(OUT_PROC_DIR, "idol_cpgs_shared.rds")
ref_centroids_file <- file.path(OUT_PROC_DIR, "ref_beta_shared_centroids.rds")
targets_file <- file.path(OUT_PROC_DIR, "targets_merged.rds")

beta_file <- file.path(OUT_PROC_DIR, "beta_matrix_sesame_batch_corrected.rds")
mval_file <- file.path(OUT_PROC_DIR, "mval_matrix_sesame_batch_corrected.rds")

stopifnot(file.exists(idol_file), file.exists(ref_centroids_file))

idol_probes <- readRDS(idol_file)
ref_centroids <- readRDS(ref_centroids_file)

# Prefer beta matrix; if missing, use mval but warn (deconv works best in beta space)
if (file.exists(beta_file)) {
  message("Loading bulk beta matrix: ", beta_file)
  bulk <- readRDS(beta_file)
  bulk_space <- "beta"
} else if (file.exists(mval_file)) {
  warning("Bulk beta matrix not found. Using mval matrix (not ideal for deconvolution): ", mval_file)
  bulk <- readRDS(mval_file)
  bulk_space <- "mval"
} else {
  stop("No bulk matrix found. Expected:\n - ", beta_file, "\n - ", mval_file)
}

# Optional targets
targets <- NULL
if (file.exists(targets_file)) {
  targets <- readRDS(targets_file)
}

# ------------------------------------------------------------
# 2) Normalize + collapse bulk IDs
# ------------------------------------------------------------
message("Normalizing SeSAMe IDs and collapsing duplicates...")
rownames(bulk) <- normalize_sesame_ids_vec(rownames(bulk))
bulk <- collapse_duplicate_rows_mean(bulk)

message("Bulk matrix after normalize+collapse: ", nrow(bulk), " CpGs x ", ncol(bulk), " samples")

# ------------------------------------------------------------
# 3) Restrict to IDOL probes and align with reference centroids
# ------------------------------------------------------------
common <- intersect(idol_probes, rownames(bulk))
message("IDOL probes present in bulk: ", length(common), " / ", length(idol_probes))

if (length(common) < 100) {
  stop("Too few IDOL probes found in bulk after normalization. ",
       "This indicates an ID mismatch or wrong files.")
}

# ref_centroids is on "shared CpGs" universe; ensure it contains these probes
if (!all(common %in% rownames(ref_centroids))) {
  missing_in_ref <- setdiff(common, rownames(ref_centroids))
  stop("Some IDOL probes are missing from reference centroids (unexpected). Missing n=",
       length(missing_in_ref))
}

Y <- bulk[common, , drop = FALSE]
R <- ref_centroids[common, , drop = FALSE]

# ------------------------------------------------------------
# 4) Deconvolve
# ------------------------------------------------------------
message("Deconvolving using ", bulk_space, " space...")
fractions <- solve_fractions_nnls(Y = Y, R = R)

# Convert to samples x celltypes for convenience
fractions_df <- as.data.frame(t(fractions))
fractions_df$sample_name <- rownames(fractions_df)

# ------------------------------------------------------------
# 5) Merge targets (if available) and save
# ------------------------------------------------------------
if (!is.null(targets) && ("sample_name" %in% names(targets))) {
  # keep only relevant cols to avoid huge merges
  keep_cols <- intersect(names(targets), c("sample_name", "patient_id", "material", "type", "source", "batch"))
  merged <- merge(
    x = targets[, keep_cols, drop = FALSE],
    y = fractions_df,
    by = "sample_name",
    all.y = TRUE,
    sort = FALSE
  )
} else {
  merged <- fractions_df
}

out_rds <- file.path(OUT_RES_DIR, "cell_fractions_idol_nnls.rds")
out_csv <- file.path(OUT_RES_DIR, "cell_fractions_idol_nnls.csv")

saveRDS(merged, out_rds)
# Round numeric columns to 5 digits
merged_csv <- merged
num_cols <- sapply(merged_csv, is.numeric)
merged_csv[, num_cols] <- round(merged_csv[, num_cols], 5)

fwrite(merged_csv, out_csv)

message("Done. Saved:")
message(" - ", out_rds)
message(" - ", out_csv)

# Also save a wide matrix (celltypes x samples) for downstream math
out_mat <- file.path(OUT_RES_DIR, "cell_fractions_idol_nnls_matrix.rds")
saveRDS(fractions, out_mat)
message(" - ", out_mat)

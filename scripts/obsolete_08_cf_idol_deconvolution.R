#!/usr/bin/env Rscript

# ============================================================
# 08_cf_idol_deconvolution.R
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
  library(nnls)
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
solve_fractions_nnls <- function(Y, R, min_probes = 50) {
  # Y: probes x samples
  # R: probes x celltypes (centroids)
  stopifnot(rownames(Y) == rownames(R))
  
  X <- matrix(NA_real_, nrow = ncol(R), ncol = ncol(Y))
  rownames(X) <- colnames(R)
  colnames(X) <- colnames(Y)
  
  Rt_full <- as.matrix(R)
  
  for (j in seq_len(ncol(Y))) {
    
    y <- Y[, j]
    keep <- is.finite(y)
    
    n_keep <- sum(keep)
    if (n_keep < min_probes) {
      # Not enough probes → leave NA (or zeros if you prefer)
      next
    }
    
    Rt <- Rt_full[keep, , drop = FALSE]
    yj <- y[keep]
    
    fit <- tryCatch(
      nnls::nnls(Rt, yj),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    x <- as.numeric(fit$x)
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
idol_file <- file.path(OUT_PROC_DIR, "idol_cpgs_shared_cg.rds")
ref_centroids_file <- file.path(OUT_PROC_DIR, "ref_beta_shared_centroids_cg.rds")

beta_file <- file.path(OUT_PROC_DIR, "cf_beta_matrix_sesame.rds")

stopifnot(file.exists(idol_file), file.exists(ref_centroids_file))

idol_probes <- readRDS(idol_file)
ref_centroids <- readRDS(ref_centroids_file)

bulk <- readRDS(beta_file)
bulk_space <- "beta"

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
message("DEBUG: R version: ", R.version.string)
message("DEBUG: data.table version: ", as.character(packageVersion("data.table")))
message("DEBUG: nnls available? ", requireNamespace("nnls", quietly = TRUE))
if (requireNamespace("nnls", quietly = TRUE)) {
  message("DEBUG: nnls version: ", as.character(packageVersion("nnls")))
} else {
  message("DEBUG: Using fallback lm.fit path (more fragile with NA / rank issues).")
}
fractions <- solve_fractions_nnls(Y = Y, R = R)

stopifnot(
  all(abs(colSums(fractions) - 1) < 1e-6 | colSums(fractions) == 0)
)

# Convert to samples x celltypes for convenience
fractions_df <- as.data.frame(t(fractions))
fractions_df$Patient <- rownames(fractions_df)
# derive Disease from Patient (ignore first 5 chars)
fractions_df$Disease <- substr(fractions_df$Patient, 6, pmin(nchar(fractions_df$Patient), 8))

out_rds <- file.path(OUT_RES_DIR, "blood_cell_fractions_idol_cellfree.rds")
out_csv <- file.path(OUT_RES_DIR, "blood_cell_fractions_idol_cellfree.csv")

saveRDS(fractions_df, out_rds)
# Round numeric columns to 5 digits
merged_csv <- fractions_df
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

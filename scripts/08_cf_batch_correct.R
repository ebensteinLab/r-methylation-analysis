#!/usr/bin/env Rscript

# ================================================================
# Script 08: Batch correction for cfDNA methylation data
#            using chromosome-wise ComBat to avoid OOM
#
# Input:
#   - results/processed/cf_mval_matrix_raw_sesame.rds
#   - results/processed/cf_targets_with_sesame.rds
#
# Output:
#   - results/processed/cf_mval_matrix_sesame_batch_corrected.rds
#   - results/processed/cf_beta_matrix_sesame_batch_corrected.rds
#   - results/processed/cf_batch_correction_metadata.rds
# ================================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(sva)
  library(sesame)
  library(sesameData)
  library(GenomicRanges)
})

options(stringsAsFactors = FALSE)
options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

OUT_PROC_DIR <- "results/processed"
dir.create(OUT_PROC_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading cfDNA SeSAMe-processed data...")

mval <- readRDS(file.path(OUT_PROC_DIR, "cf_mval_matrix_raw_sesame.rds"))
targets <- readRDS(file.path(OUT_PROC_DIR, "cf_targets_with_sesame.rds"))

message("cf M-value matrix dims: ", nrow(mval), " x ", ncol(mval))
message("Number of rows in cf targets: ", nrow(targets))

# ------------------------------------------------
# Ensure sample alignment
# ------------------------------------------------
if (!"Patient" %in% colnames(targets)) {
  stop("cf_targets_with_sesame.rds must contain a 'Patient' column")
}

targets$Patient <- as.character(targets$Patient)

mval_samples   <- colnames(mval)
target_samples <- targets$Patient

missing_in_targets <- setdiff(mval_samples, target_samples)
missing_in_mval    <- setdiff(target_samples, mval_samples)

if (length(missing_in_targets) > 0 || length(missing_in_mval) > 0) {
  
  message("---- Sample name mismatch detected ----")
  
  if (length(missing_in_targets) > 0) {
    message("Samples in mval but NOT in targets (", length(missing_in_targets), "):")
    print(missing_in_targets)
  }
  
  if (length(missing_in_mval) > 0) {
    message("Samples in targets but NOT in mval (", length(missing_in_mval), "):")
    print(missing_in_mval)
  }
  
  stop("Sample names in cf mval matrix do not match cf targets table.")
}

targets <- targets[match(colnames(mval), targets$Patient), , drop = FALSE]

stopifnot(identical(colnames(mval), targets$Patient))

# ------------------------------------------------
# Helper functions
# ------------------------------------------------
derive_disease_from_patient <- function(patient_ids) {
  substr(patient_ids, 6, pmin(nchar(patient_ids), 8))
}

pick_first_existing_column <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) return(NULL)
  hits[[1]]
}

safe_factor <- function(x, name) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "UNKNOWN"
  f <- factor(x)
  if (nlevels(f) < 2) {
    stop("Covariate '", name, "' has fewer than 2 levels after cleaning")
  }
  f
}

# ------------------------------------------------
# Define batch variable
# ------------------------------------------------
BATCH_COL_CANDIDATES <- c(
  "Batch",
  "batch",
  "Slide",
  "slide",
  "Sentrix_Id",
  "sentrix_id",
  "Chip",
  "chip",
  "Plate",
  "plate",
  "Array",
  "array",
  "Run",
  "run"
)

batch_col <- pick_first_existing_column(targets, BATCH_COL_CANDIDATES)

if (is.null(batch_col)) {
  stop(
    "Could not find a batch column in cf targets. Tried: ",
    paste(BATCH_COL_CANDIDATES, collapse = ", "),
    ". Please edit BATCH_COL_CANDIDATES or hard-code the batch column."
  )
}

batch <- safe_factor(targets[[batch_col]], batch_col)

if (anyNA(batch)) {
  stop("Batch variable contains NA values.")
}

message("Using batch column: ", batch_col)
message("Batch levels: ", paste(levels(batch), collapse = ", "))
print(sort(table(batch), decreasing = TRUE))

# ------------------------------------------------
# Build model matrix with protected biological covariates
# ------------------------------------------------
if ("Disease" %in% colnames(targets)) {
  disease_raw <- targets$Disease
  message("Using existing Disease column from cf targets")
} else {
  disease_raw <- derive_disease_from_patient(targets$Patient)
  message("Disease column not found; deriving Disease from Patient")
}

covar_df <- data.frame(
  Disease = safe_factor(disease_raw, "Disease")
)

# Optional covariates, if available and informative
for (nm in c("Sex", "Gender")) {
  if (nm %in% colnames(targets)) {
    x <- as.character(targets[[nm]])
    x[is.na(x) | trimws(x) == ""] <- "UNKNOWN"
    xf <- factor(x)
    if (nlevels(xf) >= 2) {
      covar_df[[nm]] <- xf
      message("Including categorical covariate: ", nm)
      break
    }
  }
}

if ("Age" %in% colnames(targets)) {
  age_num <- suppressWarnings(as.numeric(targets$Age))
  if (sum(!is.na(age_num)) >= 3 && sd(age_num, na.rm = TRUE) > 0) {
    covar_df$Age <- age_num
    message("Including numeric covariate: Age")
  }
}

mod <- model.matrix(~ ., data = covar_df)

message("Protected covariates in ComBat model:")
print(colnames(mod))

MAX_NA_FRAC_PROBE <- 0.20

finite_frac <- rowMeans(is.finite(mval))
keep <- finite_frac >= (1 - MAX_NA_FRAC_PROBE)

message("Keeping ", sum(keep), " / ", length(keep),
        " probes with finite fraction >= ", 1 - MAX_NA_FRAC_PROBE)

mval <- mval[keep, , drop = FALSE]

# ------------------------------------------------
# Keep only probes with finite values in all samples
# ------------------------------------------------
keep <- rowSums(is.finite(mval)) == ncol(mval)
message("Keeping ", sum(keep), " / ", length(keep), " probes with all finite values")

mval <- mval[keep, , drop = FALSE]

if (nrow(mval) < 100000) {
  warning("Fewer than 100000 probes retained after finite-value filtering")
}

# ------------------------------------------------
# Run chromosome-wise ComBat to avoid OOM
# ------------------------------------------------
if (nlevels(batch) < 2) {
  warning("Only one batch detected; skipping ComBat.")
  mval_corrected <- mval
  combat_ran <- FALSE
} else {
  message("Loading EPICv2 manifest and assigning chromosomes...")
  gr <- sesameData_getManifestGRanges("EPICv2")
  
  chr <- as.character(seqnames(gr))
  names(chr) <- names(gr)
  
  chr <- chr[rownames(mval)]
  
  if (anyNA(chr)) {
    message("Some probes are missing chromosome annotations and will be left as NA unless handled below.")
  }
  
  mval_corrected <- matrix(
    NA_real_,
    nrow = nrow(mval),
    ncol = ncol(mval),
    dimnames = dimnames(mval)
  )
  
  unique_chr <- unique(chr)
  unique_chr <- unique_chr[!is.na(unique_chr)]
  
  for (c in unique_chr) {
    message("Running ComBat on ", c)
    
    idx <- which(chr == c)
    
    if (length(idx) < 100) {
      message("Skipping ", c, " because it has fewer than 100 probes")
      next
    }
    
    X <- mval[idx, , drop = FALSE]
    
    Xc <- ComBat(
      dat = X,
      batch = batch,
      mod = mod,
      par.prior = TRUE,
      prior.plots = FALSE
    )
    
    mval_corrected[idx, ] <- Xc
    
    rm(X, Xc, idx)
    gc()
  }
  
  # Handle probes without chromosome assignment or skipped tiny chromosomes:
  not_done <- which(apply(is.na(mval_corrected), 1, all))
  if (length(not_done) > 0) {
    message("Leaving ", length(not_done), " probes unchanged because they were unassigned/skipped")
    mval_corrected[not_done, ] <- mval[not_done, , drop = FALSE]
  }
  
  combat_ran <- TRUE
  rm(gr, chr, unique_chr, not_done)
  gc()
}

# ------------------------------------------------
# Save corrected M-values
# ------------------------------------------------
mval_out <- file.path(OUT_PROC_DIR, "cf_mval_matrix_sesame_batch_corrected.rds")
saveRDS(mval_corrected, mval_out)

# ------------------------------------------------
# Convert corrected M-values back to beta values
# ------------------------------------------------
message("Converting corrected cf M-values back to betas...")
# Convert corrected M-values back to betas
# Use explicit formula to preserve matrix dimensions and dimnames
beta_corrected <- 2^mval_corrected / (2^mval_corrected + 1)

# Keep beta values in valid range
beta_corrected[beta_corrected < 0] <- 0
beta_corrected[beta_corrected > 1] <- 1

# Explicitly restore dimnames as a safety check
dimnames(beta_corrected) <- dimnames(mval_corrected)

if (!is.matrix(beta_corrected)) {
  stop("beta_corrected is not a matrix after M-value to beta conversion")
}

if (is.null(rownames(beta_corrected)) || is.null(colnames(beta_corrected))) {
  stop("beta_corrected lost rownames or colnames")
}

message("Corrected beta dims: ", nrow(beta_corrected), " x ", ncol(beta_corrected))
message("Corrected beta first rowname: ", rownames(beta_corrected)[1])
message("Corrected beta first colname: ", colnames(beta_corrected)[1])

beta_out <- file.path(OUT_PROC_DIR, "cf_beta_matrix_sesame_batch_corrected.rds")
saveRDS(beta_corrected, beta_out)

# ------------------------------------------------
# Save metadata
# ------------------------------------------------
meta_out <- file.path(OUT_PROC_DIR, "cf_batch_correction_metadata.rds")
saveRDS(
  list(
    batch_column = batch_col,
    batch_sizes = sort(table(batch), decreasing = TRUE),
    protected_covariates = colnames(mod),
    kept_probe_mask = keep,
    sample_order = colnames(mval_corrected),
    combat_ran = combat_ran
  ),
  meta_out
)

rm(mval, mval_corrected, beta_corrected, targets, covar_df, mod, batch, keep)
gc()

if (combat_ran) {
  message("cfDNA batch correction completed successfully.")
} else {
  message("cfDNA batch correction skipped (single batch).")
}

message("Saved:")
message("  ", mval_out)
message("  ", beta_out)
message("  ", meta_out)
#!/usr/bin/env Rscript

# ================================================================
# Script 04: Batch correction using ComBat (SeSAMe pipeline)
#            with biological covariates
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

message("Loading SeSAMe-processed data...")

mval <- readRDS("results/processed/mval_matrix_raw_sesame.rds")
targets <- readRDS("results/processed/targets_merged.rds")

if (is.null(dim(mval))) {
  stop("mval was loaded as a vector, not a matrix")
}

message("Number of rows in targets: ", nrow(targets))

# ------------------------------------------------
# Ensure sample alignment
# ------------------------------------------------
if (!"Patient" %in% colnames(targets)) {
  stop("targets must contain a 'Patient' column")
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
  
  stop("Sample names in mval matrix do not match targets table.")
}

targets <- targets[match(colnames(mval), targets$Patient), , drop = FALSE]

if (!identical(colnames(mval), targets$Patient)) {
  stop("Column names of mval do not match targets$Patient after alignment")
}

# ------------------------------------------------
# Helpers
# ------------------------------------------------
derive_disease_from_patient <- function(patient_ids) {
  # Keep your existing convention used elsewhere
  substr(patient_ids, 6, pmin(nchar(patient_ids), 8))
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
if (!"Sentrix_Id" %in% colnames(targets)) {
  stop("targets must contain 'Sentrix_Id' for batch correction")
}

batch <- factor(targets$Sentrix_Id)

if (anyNA(batch)) {
  stop("Batch variable contains NA values.")
}

message("Batch levels: ", paste(levels(batch), collapse = ", "))
print(sort(table(batch), decreasing = TRUE))

# ------------------------------------------------
# Keep only probes with no missing values across samples
# ------------------------------------------------
keep <- rowSums(!is.na(mval)) == length(batch)
message("Keeping ", sum(keep), " / ", length(keep), " probes with complete data before ComBat")
mval <- mval[keep, , drop = FALSE]

# ------------------------------------------------
# Build biological covariate model
# ------------------------------------------------
if ("Disease" %in% colnames(targets)) {
  disease_raw <- targets$Disease
  message("Using existing Disease column from targets")
} else {
  disease_raw <- derive_disease_from_patient(targets$Patient)
  message("Disease column not found; deriving Disease from Patient")
}

covar_df <- data.frame(
  Disease = safe_factor(disease_raw, "Disease")
)

# Optional categorical covariate
sex_col <- NULL
for (nm in c("Sex", "Gender")) {
  if (nm %in% colnames(targets)) {
    x <- as.character(targets[[nm]])
    x[is.na(x) | trimws(x) == ""] <- "UNKNOWN"
    xf <- factor(x)
    if (nlevels(xf) >= 2) {
      covar_df[[nm]] <- xf
      sex_col <- nm
      message("Including categorical covariate: ", nm)
      break
    }
  }
}

# Optional numeric covariate
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

# Optional sanity check for confounding
message("Batch x Disease table:")
print(table(batch, covar_df$Disease))

# ------------------------------------------------
# Run ComBat chromosome-wise to avoid OOM
# ------------------------------------------------
if (nlevels(batch) < 2) {
  warning("Only one batch detected; skipping ComBat.")
  mval_corrected <- mval
  combat_ran <- FALSE
} else {
  gr <- sesameData_getManifestGRanges("EPICv2")
  chr <- as.character(seqnames(gr))
  names(chr) <- names(gr)
  
  # Keep only probes in matrix
  chr <- chr[rownames(mval)]
  
  mval_corrected <- matrix(
    NA_real_,
    nrow = nrow(mval),
    ncol = ncol(mval),
    dimnames = dimnames(mval)
  )
  
  if (is.null(dim(mval_corrected))) {
    stop("mval_corrected lost dimensions immediately after initialization")
  }
  
  unique_chr <- unique(chr)
  unique_chr <- unique_chr[!is.na(unique_chr)]
  
  for (c in unique_chr) {
    message("Running ComBat on ", c)
    
    idx <- which(chr == c)
    
    # Skip tiny chromosomes
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
  
  # For unassigned / skipped probes, keep original values
  not_done <- which(apply(is.na(mval_corrected), 1, all))
  if (length(not_done) > 0) {
    message("Leaving ", length(not_done), " probes unchanged because they were unassigned/skipped")
    mval_corrected[not_done, ] <- mval[not_done, , drop = FALSE]
  }
  
  combat_ran <- TRUE
  rm(gr, chr, unique_chr, not_done)
  gc()
}

protected_covariates <- colnames(mod)

rm(mval, mod, batch, keep, covar_df)
gc()

# ------------------------------------------------
# Save output
# ------------------------------------------------
message("Saving outputs")

saveRDS(mval_corrected, "results/processed/mval_matrix_sesame_batch_corrected.rds")

# Convert corrected M-values back to beta values
message("Converting corrected M-values back to beta values")

mval_dim <- dim(mval_corrected)
mval_dimnames <- dimnames(mval_corrected)

if (is.null(mval_dim) || length(mval_dim) != 2) {
  stop("mval_corrected is not a 2D matrix before beta conversion")
}

beta_corrected <- 2^mval_corrected / (2^mval_corrected + 1)

beta_corrected <- matrix(
  beta_corrected,
  nrow = mval_dim[1],
  ncol = mval_dim[2],
  dimnames = mval_dimnames
)

beta_corrected[beta_corrected < 0] <- 0
beta_corrected[beta_corrected > 1] <- 1

message("Corrected beta dims: ", nrow(beta_corrected), " x ", ncol(beta_corrected))
message("Corrected beta first rowname: ", rownames(beta_corrected)[1])
message("Corrected beta first colname: ", colnames(beta_corrected)[1])

saveRDS(beta_corrected, "results/processed/beta_matrix_sesame_batch_corrected.rds")

# Optional metadata for reproducibility
saveRDS(
  list(
    protected_covariates = protected_covariates,
    combat_ran = combat_ran,
    n_probes_corrected = nrow(mval_corrected),
    n_samples = ncol(mval_corrected)
  ),
  "results/processed/batch_correction_metadata.rds"
)

rm(mval_corrected, beta_corrected, targets)
gc()

if (combat_ran) {
  message("Batch correction completed successfully.")
} else {
  message("Batch correction skipped (single batch).")
}


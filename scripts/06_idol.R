#!/usr/bin/env Rscript

# ============================================================
# 06_idol.R
# Train IDOL using FlowSorted.Blood.EPIC reference, restricted
# to CpGs present in your EPICv2 SeSAMe bulk matrices.
#
# Key fixes:
# - NO liftover (Illumina cg IDs are shared directly)
# - Normalize SeSAMe suffix IDs
# - Collapse duplicates after normalization
# - Restrict training universe to probes present in bulk
# ============================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(ExperimentHub)
  library(minfi)
  library(IDOL)
})

options(ExperimentHub.ask = FALSE)
set.seed(1)

OUT_PROC_DIR <- "results/processed"
dir.create(OUT_PROC_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

# cg00000109_TC21 -> cg00000109
normalize_sesame_ids_vec <- function(x) sub("_.*$", "", x)

# Collapse duplicated rownames by mean across duplicates
collapse_duplicate_rows_mean <- function(mat) {
  stopifnot(!is.null(rownames(mat)))
  rn <- rownames(mat)
  if (!anyDuplicated(rn)) return(mat)
  
  # rowsum collapses by group; works for numeric matrices
  # but we need means: sum / count
  counts <- table(rn)
  summed <- rowsum(mat, group = rn, reorder = FALSE)
  means <- summed / as.numeric(counts[rownames(summed)])
  means
}

# Load bulk matrix (beta preferred, else mval) and return normalized+collapsed rownames
load_bulk_probe_universe <- function(out_dir) {
  beta_file <- file.path(out_dir, "beta_matrix_sesame_batch_corrected.rds")
  mval_file <- file.path(out_dir, "mval_matrix_sesame_batch_corrected.rds")
  
  if (file.exists(beta_file)) {
    message("Loading bulk beta matrix: ", beta_file)
    bulk <- readRDS(beta_file)
  } else if (file.exists(mval_file)) {
    message("Loading bulk mval matrix: ", mval_file)
    bulk <- readRDS(mval_file)
  } else {
    stop("Could not find bulk beta or mval matrix in ", out_dir,
         "\nExpected one of:\n - ", beta_file, "\n - ", mval_file)
  }
  
  # Normalize and collapse duplicates
  rownames(bulk) <- normalize_sesame_ids_vec(rownames(bulk))
  bulk <- collapse_duplicate_rows_mean(bulk)
  
  message("Bulk probe universe after normalize+collapse: ", nrow(bulk))
  rownames(bulk)
}

# ------------------------------------------------------------
# 1) Load FlowSorted reference (EPIC v1)
# ------------------------------------------------------------
message("Loading FlowSorted Blood EPIC reference (EH1136)...")
eh <- ExperimentHub()
ref_rg <- eh[["EH1136"]]

ref_mset <- preprocessNoob(ref_rg)
ref_beta <- getBeta(ref_mset)
celltypes <- factor(pData(ref_mset)$CellType)

rm(ref_rg, ref_mset)
gc()

message("Reference samples: ", length(celltypes))
message("Reference cell types: ", paste(levels(celltypes), collapse = ", "))
message("Reference CpGs: ", nrow(ref_beta))

# ------------------------------------------------------------
# 2) Restrict training CpGs to those present in your EPICv2 bulk matrices
#    (by cg ID, after SeSAMe suffix normalization)
# ------------------------------------------------------------
bulk_universe <- load_bulk_probe_universe(OUT_PROC_DIR)

common <- intersect(rownames(ref_beta), bulk_universe)
message("CpGs shared (reference ∩ bulk): ", length(common))

if (length(common) < 50000) {
  stop("Too few shared CpGs (", length(common), "). ",
       "This likely indicates a probe-ID mismatch or wrong input matrices.")
}

ref_beta_shared <- ref_beta[common, , drop = FALSE]
rm(ref_beta, bulk_universe, common)
gc()

message("Training beta matrix: ", nrow(ref_beta_shared), " CpGs x ",
        ncol(ref_beta_shared), " samples")

# ------------------------------------------------------------
# 3) Save EPICv2-space reference centroids for deconvolution
#    (still cg IDs; these are what exist in your EPICv2 bulk after normalization)
# ------------------------------------------------------------
message("Computing reference centroids (shared CpGs)...")
ref_centroids <- sapply(levels(celltypes), function(ct) {
  rowMeans(ref_beta_shared[, celltypes == ct, drop = FALSE], na.rm = TRUE)
})
ref_centroids <- as.matrix(ref_centroids)

saveRDS(ref_centroids, file.path(OUT_PROC_DIR, "ref_beta_shared_centroids.rds"))
saveRDS(levels(celltypes), file.path(OUT_PROC_DIR, "ref_celltypes_levels.rds"))

rm(ref_centroids)
gc()

# ------------------------------------------------------------
# 4) Candidate DMR Finder
# ------------------------------------------------------------
message("Creating CandidateDMRFinder.v2...")
covars <- data.frame(CellType = celltypes, dummy = seq_along(celltypes))
rownames(covars) <- colnames(ref_beta_shared)

candFinder <- CandidateDMRFinder.v2(
  cellTypes       = celltypes,
  referenceBetas  = ref_beta_shared,
  referenceCovars = covars,
  M               = 150,
  equal.variance  = FALSE
)

rm(covars)
gc()

# ------------------------------------------------------------
# 5) Patch coefEsts if duplicated columns exist (common with some setups)
# ------------------------------------------------------------
if (!is.null(candFinder$coefEsts) && anyDuplicated(colnames(candFinder$coefEsts))) {
  message("Patching candFinder$coefEsts (collapsing duplicate columns by mean)...")
  coef <- candFinder$coefEsts
  labs <- colnames(coef)
  u <- unique(labs)
  coef_collapsed <- sapply(u, function(ct) rowMeans(coef[, labs == ct, drop = FALSE]))
  candFinder$coefEsts <- coef_collapsed
  rm(coef, labs, u, coef_collapsed)
  gc()
}

idol_classes <- colnames(candFinder$coefEsts)
message("IDOL classes: ", paste(idol_classes, collapse = ", "))

# ------------------------------------------------------------
# 6) Build training covariates (one-hot)
# ------------------------------------------------------------
onehot <- model.matrix(~ 0 + celltypes)
colnames(onehot) <- levels(celltypes)

trainingCovariates <- onehot[, idol_classes, drop = FALSE]
rownames(trainingCovariates) <- colnames(ref_beta_shared)

stopifnot(all(colnames(trainingCovariates) == idol_classes))

rm(onehot)
gc()

# ------------------------------------------------------------
# 7) Run IDOL optimization
# ------------------------------------------------------------
message("Running IDOLoptimize...")
idol_res <- IDOLoptimize(
  candDMRFinderObject = candFinder,
  trainingBetas       = ref_beta_shared,
  trainingCovariates  = trainingCovariates,
  libSize             = 300,
  maxIt               = 500,
  numCores            = 4
)

idol_probes <- idol_res[["IDOL Optimized Library"]]
message("IDOL selected CpGs: ", length(idol_probes))

# ------------------------------------------------------------
# 8) Save outputs
# ------------------------------------------------------------
saveRDS(idol_probes, file.path(OUT_PROC_DIR, "idol_cpgs_shared.rds"))
saveRDS(idol_res,    file.path(OUT_PROC_DIR, "idol_model_shared.rds"))
saveRDS(ref_beta_shared, file.path(OUT_PROC_DIR, "ref_beta_shared.rds"))

rm(candFinder, ref_beta_shared, trainingCovariates, idol_res, idol_probes, celltypes)
gc()

message("Done. Saved:")
message(" - results/processed/idol_cpgs_shared.rds")
message(" - results/processed/idol_model_shared.rds")
message(" - results/processed/ref_beta_shared_centroids.rds")
message(" - results/processed/ref_beta_shared.rds (large)")
message(" - results/processed/ref_celltypes_levels.rds")

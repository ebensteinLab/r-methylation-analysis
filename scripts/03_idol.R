#!/usr/bin/env Rscript

# ============================================================
# Script 03: IDOL training (genomic)
#
# Train IDOL using FlowSorted.Blood.EPIC reference,
# restricted to CpGs present in your EPICv2 SeSAMe bulk matrix.
#
# Design principles:
# - NO SeSAMe mask usage
# - NO liftover
# - Restrict by bulk probe universe only
# - EPICv2 suffix IDs handled explicitly
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
normalize_sesame_ids <- function(x) sub("_.*$", "", x)

# Collapse duplicated rownames by mean
collapse_duplicate_rows_mean <- function(mat) {
  rn <- rownames(mat)
  if (!anyDuplicated(rn)) return(mat)
  
  counts <- table(rn)
  summed <- rowsum(mat, group = rn, reorder = FALSE)
  summed / as.numeric(counts[rownames(summed)])
}

# ------------------------------------------------------------
# 1) Load bulk EPICv2 beta matrix (defines training universe)
# ------------------------------------------------------------
message("Loading bulk EPICv2 beta matrix...")
bulk_beta <- readRDS(file.path(OUT_PROC_DIR, "beta_matrix_sesame.rds"))

# Normalize suffix IDs → cg IDs
bulk_cg <- normalize_sesame_ids(rownames(bulk_beta))

# Build cg → suffix mapping (first occurrence is fine)
cg_to_full <- tapply(rownames(bulk_beta), bulk_cg, function(v) v[[1]])

# Collapse duplicated cg rows
rownames(bulk_beta) <- bulk_cg
bulk_beta <- collapse_duplicate_rows_mean(bulk_beta)

bulk_probe_universe <- rownames(bulk_beta)
message("Bulk probe universe (cg IDs): ", length(bulk_probe_universe))

rm(bulk_beta)
gc()

# ------------------------------------------------------------
# 2) Load FlowSorted Blood EPIC reference
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
# 3) Restrict reference to bulk probe universe
# ------------------------------------------------------------
common_cg <- intersect(rownames(ref_beta), bulk_probe_universe)
message("CpGs shared (reference ∩ bulk): ", length(common_cg))

if (length(common_cg) < 50000) {
  stop("Too few shared CpGs (", length(common_cg), "). Check probe IDs.")
}

ref_beta_shared <- ref_beta[common_cg, , drop = FALSE]
rm(ref_beta)
gc()

message("Training matrix: ", nrow(ref_beta_shared), " CpGs x ", ncol(ref_beta_shared), " samples")

# ------------------------------------------------------------
# 4) Save reference centroids (cg-space)
# ------------------------------------------------------------
message("Computing reference centroids...")
ref_centroids_cg <- sapply(levels(celltypes), function(ct) {
  rowMeans(ref_beta_shared[, celltypes == ct, drop = FALSE], na.rm = TRUE)
})
ref_centroids_cg <- as.matrix(ref_centroids_cg)

saveRDS(ref_centroids_cg, file.path(OUT_PROC_DIR, "ref_beta_shared_centroids_cg.rds"))
saveRDS(levels(celltypes), file.path(OUT_PROC_DIR, "ref_celltypes_levels.rds"))

# ------------------------------------------------------------
# 5) Candidate DMR Finder
# ------------------------------------------------------------
message("Creating CandidateDMRFinder.v2...")
covars <- data.frame(CellType = celltypes, dummy = seq_along(celltypes))
rownames(covars) <- colnames(ref_beta_shared)

celltypes_unique <- levels(celltypes)

candFinder <- CandidateDMRFinder.v2(
  cellTypes       = celltypes_unique,
  referenceBetas  = ref_beta_shared,
  referenceCovars = covars,
  M               = 150,
  equal.variance  = FALSE
)

rm(covars, celltypes_unique)
gc()

idol_classes <- colnames(candFinder$coefEsts)
message("IDOL classes: ", paste(idol_classes, collapse = ", "))

# ------------------------------------------------------------
# 6) Training covariates (one-hot)
# ------------------------------------------------------------
onehot <- model.matrix(~ 0 + celltypes)
colnames(onehot) <- levels(celltypes)

trainingCovariates <- onehot[, idol_classes, drop = FALSE]
rownames(trainingCovariates) <- colnames(ref_beta_shared)

rm(onehot, celltypes)
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

idol_cg <- idol_res[["IDOL Optimized Library"]]
message("IDOL selected CpGs (cg IDs): ", length(idol_cg))

# ------------------------------------------------------------
# 8) Map back to EPICv2 suffix IDs
# ------------------------------------------------------------
idol_epicv2 <- unname(cg_to_full[idol_cg])
idol_epicv2 <- idol_epicv2[!is.na(idol_epicv2)]

stopifnot(length(idol_epicv2) > 200)

message("IDOL selected CpGs (EPICv2 IDs): ", length(idol_epicv2))

# EPICv2-space centroids
centroids_epicv2 <- ref_centroids_cg
rownames(centroids_epicv2) <- unname(cg_to_full[rownames(ref_centroids_cg)])
centroids_epicv2 <- centroids_epicv2[!is.na(rownames(centroids_epicv2)), , drop = FALSE]

# ------------------------------------------------------------
# 9) Save outputs
# ------------------------------------------------------------
saveRDS(idol_cg, file.path(OUT_PROC_DIR, "idol_cpgs_shared_cg.rds"))
saveRDS(idol_epicv2, file.path(OUT_PROC_DIR, "idol_cpgs_epicv2.rds"))
saveRDS(centroids_epicv2, file.path(OUT_PROC_DIR, "ref_beta_epicv2_centroids.rds"))
saveRDS(idol_res, file.path(OUT_PROC_DIR, "idol_model_shared.rds"))

rm(ref_beta_shared, candFinder, trainingCovariates, idol_res)
gc()

message("IDOL training completed successfully.")

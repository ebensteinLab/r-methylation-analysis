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
options(AnnotationHub.ask = FALSE)

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

# Collapse duplicated cg rows
rownames(bulk_beta) <- bulk_cg
bulk_beta <- collapse_duplicate_rows_mean(bulk_beta)

bulk_probe_universe <- rownames(bulk_beta)
message("Bulk probe universe (cg IDs): ", length(bulk_probe_universe))

rm(bulk_beta, bulk_cg)
gc()

# ------------------------------------------------------------
# 2) Load FlowSorted Blood EPIC reference
# ------------------------------------------------------------
message("Loading FlowSorted BloodExtended EPIC reference (EH5425)...")
eh <- ExperimentHub()
ref_rg <- eh[["EH5425"]]

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
rm(bulk_probe_universe)

if (length(common_cg) < 50000) {
  stop("Too few shared CpGs (", length(common_cg), "). Check probe IDs.")
}

ref_beta_shared <- ref_beta[common_cg, , drop = FALSE]
rm(ref_beta, common_cg)
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

saveRDS(ref_centroids_cg, file.path(OUT_PROC_DIR, "ref_beta_shared_centroids_cg_12.rds"))
saveRDS(levels(celltypes), file.path(OUT_PROC_DIR, "ref_celltypes_levels_12.rds"))

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

k_values <- seq(100, 200, by = 100)
seeds <- 1:3
rmse_results <- data.frame(K = integer(), RMSE = numeric())
idol_results_list <- list()

message("Running IDOLoptimize...")
for (k in k_values) {
  for (s in seeds) {
    seed_val <- sample.int(.Machine$integer.max, 1)
    set.seed(seed_val)
    message("Running IDOL with K = ", k, " and seed ", seed_val, " ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    
    idol_res <- IDOLoptimize(
      candDMRFinderObject = candFinder,
      trainingBetas       = ref_beta_shared,
      trainingCovariates  = trainingCovariates,
      libSize             = k,
      maxIt               = 500,
      numCores            = 4
    )
    
    # Extract RMSE
    rmse <- NA
    
    if ("RMSE" %in% names(idol_res)) {
      rmse <- idol_res$RMSE
    } else if ("OptResults" %in% names(idol_res)) {
      rmse <- min(idol_res$OptResults$RMSE, na.rm = TRUE)
    }
    
    # Store RMSE with seed
    rmse_results <- rbind(rmse_results, data.frame(K = k, seed = seed_val, RMSE = rmse))
    
    # Store full model per (K, seed)
    if (is.null(idol_results_list[[as.character(k)]])) {
      idol_results_list[[as.character(k)]] <- list()
    }
    
    idol_results_list[[as.character(k)]][[as.character(s)]] <- idol_res
    
    current_min <- min(rmse_results$RMSE, na.rm = TRUE)
    
    message("Current global min RMSE: ", signif(current_min, 5))
    
    rm(idol_res)
  }
}

rm(ref_beta_shared, ref_centroids_cg)

saveRDS(rmse_results, file.path(OUT_PROC_DIR, "idol_k_sweep_rmse.rds"))
saveRDS(idol_results_list, file.path(OUT_PROC_DIR, "idol_k_sweep_models.rds"))
write.csv(rmse_results, file = "results/idol/rmse_results_12.csv")

rm(rmse_results, idol_results_list)
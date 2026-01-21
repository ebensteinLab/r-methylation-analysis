#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(sesameData)
  library(sesame)
})

sesameDataCache()

# Prevent any hub prompts 
options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

message("Loading merged sample sheet...")
targets <- readRDS("results/processed/targets_merged.rds")

dir.create("results/processed", recursive = TRUE, showWarnings = FALSE)

process_sample <- function(i, basename) {
  
  message("Processing ", i, ": ", basename)
  
  # 1. Read IDATs
  sdf <- sesame::readIDATpair(basename, verbose = TRUE)

  # 2. Full SeSAMe preprocessing (NOOB + dye bias + probe QC)
  sdf <- sesame::prepSesame(sdf)
  betas <- getBetas(sdf)
  
  # 3. Apply Capper-style probe filtering at the SigDF level
  qc <- sesame::sesameQC_calcStats(sdf)
  qc_df <- sesame::sesameQCtoDF(qc)

  # qc_df has one row per probe
  # good == TRUE means probe passes all SeSAMe filters
  mask <- qc_df$good
  names(mask) <- qc_df$probe_id
  mask <- mask & grepl("^cg", names(mask))
  mval_raw <- BetaValueToMValue(betas)
  
  betas_masked <- betas
  betas_masked[!mask] <- NA
  mval_masked <- BetaValueToMValue(betas_masked)
  
  rm(sdf, betas, qc, qc_df)
  gc()
  list(beta = betas_masked, mval = mval_masked, mval_raw = mval_raw, mask = mask)
}

message("Starting preprocessing of all samples...")

message("Number of rows in targets: ", nrow(targets))

basenames  <- targets$Basename
sample_ids <- make.unique(as.character(targets$sample_name))

beta_list <- setNames(vector("list", length(basenames)), sample_ids)
mval_list_raw <- setNames(vector("list", length(basenames)), sample_ids)
mval_list_masked <- setNames(vector("list", length(basenames)), sample_ids)
mask_list <- setNames(vector("list", length(basenames)), sample_ids)

for (i in seq_along(basenames)) {
  res <- process_sample(i, basenames[i])
  beta_list[[sample_ids[i]]]        <- res$beta
  mval_list_raw[[sample_ids[i]]]    <- res$mval_raw
  mval_list_masked[[sample_ids[i]]] <- res$mval
  mask_list[[sample_ids[i]]] <- res$mask
}

# ---- Validation ----
if (any(!vapply(beta_list, is.numeric, logical(1)))) {
  bad <- names(beta_list)[!vapply(beta_list, is.numeric, logical(1))]
  stop("Non-numeric beta vectors: ", paste(bad, collapse = ", "))
}

all_masks <- Reduce("&", mask_list)
rm(mask_list)
saveRDS(all_masks, "results/processed/final_probe_mask.rds")
rm(all_masks)
gc()

message("Combining results into matrices...")

all_probes <- Reduce(intersect, lapply(beta_list, names))
if (length(all_probes) == 0) stop("No common probes found across samples")

message("Saving outputs...")

beta_matrix <- do.call(cbind, lapply(beta_list, function(x) x[all_probes]))
rm(beta_list)
gc()
colnames(beta_matrix) <- sample_ids
rownames(beta_matrix) <- all_probes
saveRDS(beta_matrix, "results/processed/beta_matrix_sesame.rds")
rm(beta_matrix)
gc()

mval_matrix_raw <- do.call(cbind, lapply(mval_list_raw, function(x) x[all_probes]))
rm(mval_list_raw)
gc()
colnames(mval_matrix_raw) <- sample_ids
rownames(mval_matrix_raw) <- all_probes
saveRDS(mval_matrix_raw, "results/processed/mval_matrix_raw_sesame.rds")
rm(mval_matrix_raw)
gc()

mval_matrix_masked <- do.call(cbind, lapply(mval_list_masked, function(x) x[all_probes]))
rm(mval_list_masked)
gc()
colnames(mval_matrix_masked) <- sample_ids
rownames(mval_matrix_masked) <- all_probes
saveRDS(mval_matrix_masked, "results/processed/mval_matrix_sesame.rds")
rm(mval_matrix_masked)
gc()


message("Completed Script 02 using SeSAMe.")

#!/usr/bin/env Rscript

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

process_sample <- function(basename) {
  
  message("Processing: ", basename)
  
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
  betas[!mask] <- NA
  
  mval  <- sesame::BetaValueToMValue(betas)
  
  list(beta = betas, mval = mval)
}

message("Starting preprocessing of all samples...")

basenames  <- targets$Basename
sample_ids <- make.unique(as.character(targets$sample_name))

beta_list <- setNames(vector("list", length(basenames)), sample_ids)
mval_list <- setNames(vector("list", length(basenames)), sample_ids)

for (i in seq_along(basenames)) {
  res <- process_sample(basenames[i])
  beta_list[[sample_ids[i]]] <- res$beta
  mval_list[[sample_ids[i]]] <- res$mval
}

# ---- Validation ----
if (any(!vapply(beta_list, is.numeric, logical(1)))) {
  bad <- names(beta_list)[!vapply(beta_list, is.numeric, logical(1))]
  stop("Non-numeric beta vectors: ", paste(bad, collapse = ", "))
}

message("Combining results into matrices...")

all_probes <- Reduce(intersect, lapply(beta_list, names))
if (length(all_probes) == 0) stop("No common probes found across samples")

beta_matrix <- do.call(cbind, lapply(beta_list, function(x) x[all_probes]))
mval_matrix <- do.call(cbind, lapply(mval_list, function(x) x[all_probes]))

colnames(beta_matrix) <- sample_ids
colnames(mval_matrix) <- sample_ids
rownames(beta_matrix) <- all_probes
rownames(mval_matrix) <- all_probes

message("Saving outputs...")
saveRDS(beta_matrix, "results/processed/beta_matrix_sesame.rds")
saveRDS(mval_matrix, "results/processed/mval_matrix_sesame.rds")
saveRDS(targets,      "results/processed/targets_with_sesame.rds")

message("Completed Script 02 using SeSAMe.")

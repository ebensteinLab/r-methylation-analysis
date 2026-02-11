#!/usr/bin/env Rscript

# ================================================================
# Script 07: SeSAMe preprocessing in batches + final merge
# ================================================================

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

suppressPackageStartupMessages({
  library(sesame)
  library(sesameData)
  library(qs)
})

sesameDataCache()

options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

message("Loading sample metadata...")
cf_targets <- readRDS("results/processed/cf_targets_merged.rds")

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
BATCH_SIZE <- 100
OUT_BATCH_DIR <- "results/batches"
OUT_PROC_DIR  <- "results/processed"

dir.create(OUT_BATCH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PROC_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Split into batches
# ------------------------------------------------------------
n <- nrow(cf_targets)
batch_index <- ceiling(seq_len(n) / BATCH_SIZE)
batches <- split(cf_targets, batch_index)

saveRDS(cf_targets, file.path(OUT_PROC_DIR, "cf_targets_with_sesame.rds"))
rm(cf_targets, batch_index)
gc()

message("Total samples: ", n)
message("Batch size: ", BATCH_SIZE)
message("Total batches: ", length(batches))

# ------------------------------------------------------------
# Function to process one batch
# ------------------------------------------------------------
process_batch <- function(batch_cf_targets, batch_id) {
  
  message("\n========== Processing batch ", batch_id, " ==========")
  
  beta_list <- list()
  mval_raw_list <- list()

  for (i in seq_len(nrow(batch_cf_targets))) {
    
    basename  <- batch_cf_targets$Basename[i]
    sample_id <- batch_cf_targets$Patient[i]
    
    message("[Batch ", batch_id, "] Sample ", i, "/", nrow(batch_cf_targets),
            ": ", sample_id)
    
    sdf <- sesame::readIDATpair(basename, verbose = FALSE)
    sdf <- sesame::prepSesame(sdf)
    
    betas <- sesame::getBetas(sdf)
    mval_raw <- sesame::BetaValueToMValue(betas)
    
    beta_list[[sample_id]] <- betas
    mval_raw_list[[sample_id]] <- mval_raw

    rm(sdf, betas, mval_raw)
    gc(verbose = FALSE)
  }
  
  message("Saving batch ", batch_id)
  
  qsave(beta_list, file.path(OUT_BATCH_DIR, sprintf("cf_beta_batch_%03d.qs", batch_id)), preset="balanced")
  qsave(mval_raw_list, file.path(OUT_BATCH_DIR, sprintf("cf_mval_raw_batch_%03d.qs", batch_id)), preset="balanced")

  rm(beta_list, mval_raw_list)
  gc(verbose = FALSE)
  
  message("Batch ", batch_id, " saved.")
}


# ------------------------------------------------------------
# Run all batches (skip already processed)
# ------------------------------------------------------------
for (b in seq_along(batches)) {
  beta_file <- file.path(OUT_BATCH_DIR, sprintf("cf_beta_batch_%03d.qs", b))
  mval_raw_file <- file.path(OUT_BATCH_DIR, sprintf("cf_mval_raw_batch_%03d.qs", b))

  if (file.exists(beta_file) && file.exists(mval_raw_file)) {
    message("Skipping batch ", b, " (already exists)")
    next
  }
  
  process_batch(batches[[b]], b)
}

rm(batches, beta_file, mval_raw_file)
gc()

# ============================================================
# MERGE ALL BATCHES
# ============================================================

message("\n========== Merging batches ==========")

beta_files     <- list.files(OUT_BATCH_DIR, pattern="^cf_beta_batch_.*\\.qs$", full.names=TRUE)
mval_raw_files <- list.files(OUT_BATCH_DIR, pattern="^cf_mval_raw_batch_.*\\.qs$", full.names=TRUE)

# ---- Find common probes per batch ----
probe_sets <- lapply(beta_files, function(f) {
  batch <- qread(f)
  Reduce(intersect, lapply(batch, names))
})

all_probes <- Reduce(intersect, probe_sets)
stopifnot(length(all_probes) > 100000)
saveRDS(all_probes, file.path(OUT_PROC_DIR, "cf_common_probes.rds"))
message("Common probes: ", length(all_probes))

rm(probe_sets)
gc()

# ---- Build BETA matrix incrementally ----
sample_ids <- c()
beta_matrix <- matrix(NA, nrow = length(all_probes), ncol = 0,
                      dimnames = list(all_probes, NULL))

for (f in beta_files) {
  batch <- qread(f)
  for (sid in names(batch)) {
    sample_ids <- c(sample_ids, sid)
    vec <- batch[[sid]][all_probes]
    beta_matrix <- cbind(beta_matrix, vec)
  }
  rm(batch)
  gc()
}

colnames(beta_matrix) <- sample_ids
saveRDS(beta_matrix, file.path(OUT_PROC_DIR, "cf_beta_matrix_sesame.rds"))
rm(beta_matrix)
gc()

# ---- Build RAW M-value matrix incrementally ----
sample_ids <- c()
mval_raw_matrix <- matrix(NA, nrow=length(all_probes), ncol=0,
                          dimnames=list(all_probes, NULL))

for (f in mval_raw_files) {
  batch <- qread(f)
  for (sid in names(batch)) {
    sample_ids <- c(sample_ids, sid)
    vec <- batch[[sid]][all_probes]
    mval_raw_matrix <- cbind(mval_raw_matrix, vec)
  }
  rm(batch, vec)
  gc()
}

colnames(mval_raw_matrix) <- sample_ids
saveRDS(mval_raw_matrix, file.path(OUT_PROC_DIR, "cf_mval_matrix_raw_sesame.rds"))
rm(mval_raw_matrix, all_probes)
gc()

message("Preprocessing completed successfully.")

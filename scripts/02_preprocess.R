#!/usr/bin/env Rscript

# ================================================================
# Script 02: SeSAMe preprocessing in batches + final merge
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
targets <- readRDS("results/processed/targets_merged.rds")

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
n <- nrow(targets)
batch_index <- ceiling(seq_len(n) / BATCH_SIZE)
batches <- split(targets, batch_index)

saveRDS(targets, file.path(OUT_PROC_DIR, "targets_with_sesame.rds"))
rm(targets, batch_index)
gc()

message("Total samples: ", n)
message("Batch size: ", BATCH_SIZE)
message("Total batches: ", length(batches))

# ------------------------------------------------------------
# Function to process one batch
# ------------------------------------------------------------
process_batch <- function(batch_targets, batch_id) {
  
  message("\n========== Processing batch ", batch_id, " ==========")
  
  beta_list <- list()
  mval_list <- list()
  mval_raw_list <- list()
  mask_list <- list()
  
  for (i in seq_len(nrow(batch_targets))) {
    
    basename  <- batch_targets$Basename[i]
    sample_id <- batch_targets$sample_name[i]
    
    message("[Batch ", batch_id, "] Sample ", i, "/", nrow(batch_targets),
            ": ", sample_id)
    
    sdf <- sesame::readIDATpair(basename, verbose = FALSE)
    sdf <- sesame::prepSesame(sdf)
    
    betas <- sesame::getBetas(sdf)
    mval_raw <- sesame::BetaValueToMValue(betas)
    
    qc <- sesame::sesameQC_calcStats(sdf)
    qc_df <- sesame::sesameQCtoDF(qc)
    
    mask <- qc_df$good & grepl("^cg", qc_df$probe_id)
    names(mask) <- qc_df$probe_id
    
    betas_masked <- betas
    betas_masked[!mask] <- NA
    mval_masked <- sesame::BetaValueToMValue(betas_masked)
    
    beta_list[[sample_id]] <- betas_masked
    mval_list[[sample_id]] <- mval_masked
    mval_raw_list[[sample_id]] <- mval_raw
    mask_list[[sample_id]] <- mask
    
    rm(sdf, betas, qc, qc_df, mask, betas_masked, mval_masked, mval_raw)
    gc(verbose = FALSE)
  }
  
  message("Saving batch ", batch_id)
  
  qsave(beta_list, file.path(OUT_BATCH_DIR, sprintf("g_beta_batch_%03d.qs", batch_id)), preset="balanced")
  qsave(mval_list, file.path(OUT_BATCH_DIR, sprintf("g_mval_batch_%03d.qs", batch_id)), preset="balanced")
  qsave(mval_raw_list, file.path(OUT_BATCH_DIR, sprintf("g_mval_raw_batch_%03d.qs", batch_id)), preset="balanced")
  qsave(mask_list, file.path(OUT_BATCH_DIR, sprintf("g_mask_batch_%03d.qs", batch_id)), preset="balanced")
  
  rm(beta_list, mval_list, mval_raw_list, mask_list)
  gc(verbose = FALSE)
  
  message("Batch ", batch_id, " saved.")
}


# ------------------------------------------------------------
# Run all batches (skip already processed)
# ------------------------------------------------------------
for (b in seq_along(batches)) {
  beta_file <- file.path(OUT_BATCH_DIR, sprintf("g_beta_batch_%03d.qs", b))
  mval_file <- file.path(OUT_BATCH_DIR, sprintf("g_mval_batch_%03d.qs", b))
  mval_raw_file <- file.path(OUT_BATCH_DIR, sprintf("g_mval_raw_batch_%03d.qs", b))
  mask_file <- file.path(OUT_BATCH_DIR, sprintf("g_mask_batch_%03d.qs", b))
  
  if (file.exists(beta_file) && file.exists(mval_file)) {
    message("Skipping batch ", b, " (already exists)")
    next
  }
  
  process_batch(batches[[b]], b)
}

rm(batches, beta_file, mval_file)
gc()

# ============================================================
# MERGE ALL BATCHES
# ============================================================

message("\n========== Merging batches ==========")

# Load batches one by one
beta_files <- list.files(OUT_BATCH_DIR, pattern = "^g_beta_batch_.*qs$", full.names = TRUE)
beta_all <- list()
for (f in beta_files) {
  message("Loading ", f)
  beta_all <- c(beta_all, qread(f))
}
rm(beta_files)
gc()

# Find common probes
message("Finding common probes...")
all_probes <- Reduce(intersect, lapply(beta_all, names))
stopifnot(length(all_probes) > 100000)

message("Common probes: ", length(all_probes))

# Build matrices
message("Building beta matrix...")
beta_matrix <- do.call(cbind, lapply(beta_all, function(x) x[all_probes]))
sample_ids <- names(beta_all)
rm(beta_all)
gc()

colnames(beta_matrix) <- sample_ids
rownames(beta_matrix) <- all_probes
saveRDS(beta_matrix, file.path(OUT_PROC_DIR, "beta_matrix_sesame.rds"))
rm(beta_matrix)
gc()

message("Building RAW M-value matrix...")

mval_raw_files <- list.files(OUT_BATCH_DIR, pattern = "^g_mval_raw_batch_.*qs$", full.names = TRUE)
mval_raw_all <- list()

for (f in mval_raw_files) {
  message("Loading ", f)
  mval_raw_all <- c(mval_raw_all, qread(f))
}

rm(mval_raw_files)
gc()

mval_raw_matrix <- do.call(cbind, lapply(mval_raw_all, function(x) x[all_probes]))
rm(mval_raw_all)
gc()

colnames(mval_raw_matrix) <- sample_ids
rownames(mval_raw_matrix) <- all_probes
saveRDS(mval_raw_matrix, file.path(OUT_PROC_DIR, "mval_matrix_raw_sesame.rds"))

rm(mval_raw_matrix)
gc()

message("Building M-value matrix...")
mval_files <- list.files(OUT_BATCH_DIR, pattern = "^g_mval_batch_.*qs$", full.names = TRUE)
mval_all <- list()
for (f in mval_files) {
  message("Loading ", f)
  mval_all <- c(mval_all, qread(f))
}
rm(mval_files)
gc()

mval_matrix <- do.call(cbind, lapply(mval_all, function(x) x[all_probes]))
rm(mval_all)
gc()

colnames(mval_matrix) <- sample_ids
rownames(mval_matrix) <- all_probes
saveRDS(mval_matrix, file.path(OUT_PROC_DIR, "mval_matrix_sesame.rds"))
rm(mval_matrix, sample_ids)
gc()

message("Building mask matrix...")

mask_files <- list.files(OUT_BATCH_DIR, pattern = "^g_mask_batch_.*qs$", full.names = TRUE)

mask_all <- list()
for (f in mask_files) {
  message("Loading ", f)
  mask_all <- c(mask_all, qread(f))
}
rm(mask_files)
gc()

# Build matrix
mask_matrix <- do.call(cbind, lapply(mask_all, function(x) {
  v <- x[all_probes]
  v[is.na(v)] <- FALSE
  as.integer(v)
}))

colnames(mask_matrix) <- names(mask_all)
rownames(mask_matrix) <- all_probes
saveRDS(mask_matrix, file.path(OUT_PROC_DIR, "mask_matrix_sesame.rds"))

rm(mask_all, mask_matrix, all_probes)
gc()

message("Preprocessing completed successfully.")

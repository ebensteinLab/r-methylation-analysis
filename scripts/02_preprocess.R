#!/usr/bin/env Rscript

hub_dir <- path.expand("~/.cache/R/ExperimentHub")
dir.create(hub_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(EXPERIMENT_HUB_CACHE = hub_dir)

suppressPackageStartupMessages({
  library(sesameData)
  sesameDataCache()
  library(sesame)
})

# Optional: avoids any interactive prompts in some Bioc hubs
options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

message("Loading merged sample sheet...")
targets <- readRDS("results/processed/targets_merged.rds")

# Ensure output dir exists
dir.create("results/processed", recursive = TRUE, showWarnings = FALSE)

process_sample <- function(basename) {
  
  message("Processing sample: ", basename)
  
  betas <- openSesame(
    basename,
    mask = TRUE,
    prep_args = list(mask = c("quality", "snp5", "sex"))
  )
  
  if (!is.numeric(betas) || is.null(names(betas))) {
    stop("openSesame() failed for ", basename)
  }
  
  mval <- BetaValueToMValue(betas)
  
  list(
    beta = betas,
    mval = mval
  )
}

message("Starting preprocessing of all samples...")

basenames  <- targets$Basename
sample_ids <- make.unique(as.character(targets$sample_name))

beta_list <- vector("list", length(basenames))
mval_list <- vector("list", length(basenames))
names(beta_list) <- sample_ids
names(mval_list) <- sample_ids

for (i in seq_along(basenames)) {
  res <- process_sample(basenames[i])
  beta_list[[sample_ids[i]]] <- res$beta
  mval_list[[sample_ids[i]]] <- res$mval
}

# Validate
if (any(!vapply(beta_list, is.numeric, logical(1)))) {
  bad <- names(beta_list)[!vapply(beta_list, is.numeric, logical(1))]
  stop("Some beta_list entries are not numeric: ", paste(bad, collapse = ", "))
}

message("Combining results into matrices...")
all_probes <- Reduce(intersect, lapply(beta_list, names))
if (length(all_probes) == 0) stop("No common probes found across samples.")

beta_matrix <- do.call(cbind, lapply(beta_list, function(x) x[all_probes]))
mval_matrix <- do.call(cbind, lapply(mval_list, function(x) x[all_probes]))

colnames(beta_matrix) <- sample_ids
colnames(mval_matrix) <- sample_ids
rownames(beta_matrix) <- all_probes
rownames(mval_matrix) <- all_probes

message("Saving RDS outputs...")
saveRDS(beta_matrix, "results/processed/beta_matrix_sesame.rds")
saveRDS(mval_matrix, "results/processed/mval_matrix_sesame.rds")
saveRDS(targets,      "results/processed/targets_with_sesame.rds")

message("Completed Script 02 using SeSAMe.")

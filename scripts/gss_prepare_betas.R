#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sesame)
  library(data.table)
})

# ------------------------------------------------
# Inputs
# ------------------------------------------------
META_FILE <- "gss/mixtures_metadata.csv"

IDAT_DIRS <- c(
  "gss/idat/GSE182379",
  "gss/idat/GSE167998"
)

OUT_FILE <- "results/deconvolution/comparison/memmix_beta_for_python.csv"

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helpers
# ------------------------------------------------
extract_gsm <- function(x) {
  sub("^([^_]+)_.*$", "\\1", basename(x))
}

collapse_beta_to_cg <- function(beta_mat) {
  cg_ids <- sub("_.*$", "", rownames(beta_mat))
  
  if (!anyDuplicated(cg_ids)) {
    rownames(beta_mat) <- cg_ids
    return(beta_mat)
  }
  
  message("Collapsing duplicated EPICv2 CpGs to cg-only IDs by row mean...")
  
  beta_dt <- as.data.table(beta_mat, keep.rownames = "probe")
  beta_dt[, cg := sub("_.*$", "", probe)]
  beta_dt[, probe := NULL]
  
  collapsed <- beta_dt[, lapply(.SD, mean, na.rm = TRUE), by = cg]
  
  mat <- as.matrix(collapsed[, -"cg"])
  rownames(mat) <- collapsed$cg
  mode(mat) <- "numeric"
  
  mat
}

# ------------------------------------------------
# Find IDAT files
# ------------------------------------------------
message("Finding IDAT files...")

idat_files <- unlist(lapply(IDAT_DIRS, function(d) {
  list.files(d, pattern = "_Red.idat.gz$", full.names = TRUE)
}))

basenames_full <- sub("_Red.idat.gz$", "", idat_files)

idat_df <- data.frame(
  basename = basenames_full,
  gsm = extract_gsm(basenames_full),
  stringsAsFactors = FALSE
)

message("Found ", nrow(idat_df), " IDAT pairs")

# ------------------------------------------------
# Load metadata and match to IDATs
# ------------------------------------------------
meta <- fread(META_FILE)

if (!"gsm" %in% colnames(meta)) {
  stop("Metadata file must contain column: gsm")
}

if (!"Sample_title" %in% colnames(meta)) {
  stop("Metadata file must contain column: Sample_title")
}

meta_merged <- merge(meta, idat_df, by = "gsm")

message("Matched ", nrow(meta_merged), " metadata rows to IDATs")

if (nrow(meta_merged) == 0) {
  stop("No metadata rows matched IDAT files")
}

# Optional: keep only mixture samples if metadata contains others
meta_merged <- meta_merged[grepl("^memmix", Sample_title)]

message("Keeping ", nrow(meta_merged), " memmix samples")

# ------------------------------------------------
# Read IDATs and build beta matrix
# ------------------------------------------------
beta_list <- list()

for (i in seq_len(nrow(meta_merged))) {
  message(sprintf("Processing %d / %d: %s",
                  i, nrow(meta_merged), meta_merged$Sample_title[i]))
  
  s <- readIDATpair(meta_merged$basename[i])
  s <- prepSesame(s, "QCDPB")
  beta <- getBetas(s)
  
  beta_list[[meta_merged$Sample_title[i]]] <- beta
}

common_probes <- Reduce(intersect, lapply(beta_list, names))

message("Common probes across samples: ", length(common_probes))

beta_mat <- do.call(cbind, lapply(beta_list, function(x) x[common_probes]))
colnames(beta_mat) <- names(beta_list)
rownames(beta_mat) <- common_probes

message("Raw beta matrix dims: ", paste(dim(beta_mat), collapse = " x "))

# ------------------------------------------------
# Collapse EPICv2 probe IDs to cg-only IDs
# ------------------------------------------------
beta_cg <- collapse_beta_to_cg(beta_mat)

message("cg-level beta matrix dims: ", paste(dim(beta_cg), collapse = " x "))

# ------------------------------------------------
# Remove rows with all NA
# ------------------------------------------------
keep <- rowSums(!is.na(beta_cg)) > 0
beta_cg <- beta_cg[keep, , drop = FALSE]

message("After removing all-NA rows: ", paste(dim(beta_cg), collapse = " x "))

# ------------------------------------------------
# Write CSV for Python deconvolve.py
# Required format:
# acc,memmix13,memmix14,...
# cg...,0.1,0.2,...
# ------------------------------------------------
out_df <- data.frame(
  acc = rownames(beta_cg),
  beta_cg,
  check.names = FALSE
)

write.csv(
  out_df,
  OUT_FILE,
  row.names = FALSE,
  quote = FALSE
)

message("Wrote: ", OUT_FILE)

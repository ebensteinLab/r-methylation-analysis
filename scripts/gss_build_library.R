#!/usr/bin/env Rscript

# ================================================================
# Build EPIC beta matrix from GEO Salas IDATs for GSS marker selection
# ================================================================

suppressPackageStartupMessages({
  library(sesame)
  library(sesameData)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
})

options(stringsAsFactors = FALSE)
options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

needed_sesame_resources <- c(
  "EPIC.address",
  "KYCG.EPIC.Mask.20220123"
)

for (res in needed_sesame_resources) {
  message("Caching sesame resource: ", res)
  try(sesameDataCache(res), silent = FALSE)
}

# ------------------------------------------------
# User-configurable paths
# ------------------------------------------------
IDAT_ROOT <- "gss/idat"
META_CSV  <- "gss/reference_like_samples_only.csv"
OUT_DIR   <- "gss/processed"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helpers
# ------------------------------------------------

message2 <- function(...) {
  cat(sprintf(...), "\n")
}

normalize_sample_name <- function(x) {
  x |>
    basename() |>
    str_remove("\\.idat(\\.gz)?$") |>
    str_remove("_(Grn|Red)$")
}

extract_basename_from_idat <- function(path) {
  # Return full basename without color suffix / extension
  # Example:
  # /a/b/1234567890_R01C01_Grn.idat.gz -> /a/b/1234567890_R01C01
  fn <- basename(path)
  fn <- str_remove(fn, "\\.idat(\\.gz)?$")
  fn <- str_remove(fn, "_(Grn|Red)$")
  file.path(dirname(path), fn)
}

find_idat_pairs <- function(idat_root) {
  files <- list.files(idat_root, pattern = "\\.idat(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    stop("No IDAT files found under: ", idat_root)
  }
  
  bases <- unique(vapply(files, extract_basename_from_idat, character(1)))
  df <- data.frame(
    basename = bases,
    grn = paste0(bases, "_Grn.idat"),
    grn_gz = paste0(bases, "_Grn.idat.gz"),
    red = paste0(bases, "_Red.idat"),
    red_gz = paste0(bases, "_Red.idat.gz"),
    stringsAsFactors = FALSE
  )
  
  df$grn_exists <- file.exists(df$grn) | file.exists(df$grn_gz)
  df$red_exists <- file.exists(df$red) | file.exists(df$red_gz)
  
  df <- df[df$grn_exists & df$red_exists, , drop = FALSE]
  
  if (nrow(df) == 0) {
    stop("No complete Grn/Red IDAT pairs found.")
  }
  
  df$sample_stub <- normalize_sample_name(df$basename)
  df
}

read_metadata <- function(meta_csv) {
  meta <- fread(meta_csv, data.table = FALSE)
  
  required_cols <- c("gsm", "gse")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Metadata file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Make sure expected optional cols exist
  for (cc in c("cell_type", "sample_title", "source_name", "characteristics_text", "all_text")) {
    if (!cc %in% colnames(meta)) meta[[cc]] <- NA_character_
  }
  
  meta
}

extract_gsm_from_text <- function(x) {
  m <- str_extract(x, "GSM\\d+")
  ifelse(is.na(m), "", m)
}

match_idats_to_metadata <- function(idat_df, meta_df) {
  # First try direct GSM match from basename/path
  idat_df$path_text <- paste(idat_df$basename, idat_df$sample_stub)
  idat_df$gsm_from_path <- extract_gsm_from_text(idat_df$path_text)
  
  matched <- idat_df |>
    left_join(meta_df, by = c("gsm_from_path" = "gsm"))
  
  # If direct match fails, try title/source/characteristics text match against stub
  need_match <- is.na(matched$gse)
  
  if (any(need_match)) {
    meta_df2 <- meta_df |>
      mutate(
        meta_text = paste(
          ifelse(is.na(sample_title), "", sample_title),
          ifelse(is.na(source_name), "", source_name),
          ifelse(is.na(characteristics_text), "", characteristics_text),
          ifelse(is.na(all_text), "", all_text),
          sep = " | "
        )
      )
    
    for (i in which(need_match)) {
      stub <- matched$sample_stub[i]
      
      # Escape regex-sensitive chars
      stub_re <- str_replace_all(stub, "([\\.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1")
      
      hits <- which(str_detect(meta_df2$meta_text, fixed(stub)) |
                      str_detect(meta_df2$meta_text, regex(stub_re, ignore_case = TRUE)))
      
      if (length(hits) == 1) {
        matched[i, colnames(meta_df)] <- meta_df[hits, colnames(meta_df)]
      }
    }
  }
  
  matched$matched <- !is.na(matched$gse)
  matched
}

prep_sesame_one <- function(basename_no_suffix, prep = "QCDPB") {
  # Read raw IDATs as a SigDF
  s <- readIDATpair(basename_no_suffix)
  
  # Apply SeSAMe preprocessing
  s <- prepSesame(s, prep = prep)
  
  # Extract beta and M values
  betas <- getBetas(s)
  mvals <- BetaValueToMValue(betas)
  
  list(beta = betas, mval = mvals)
}

collapse_duplicate_probes <- function(mat, fun = c("mean", "median")) {
  fun <- match.arg(fun)
  probe_ids <- rownames(mat)
  
  # Normalize EPICv2 replicate suffixes if present:
  # cg00000109_TC21 -> cg00000109
  norm_ids <- str_extract(probe_ids, "^cg\\d+")
  
  # Keep original if pattern not found
  norm_ids[is.na(norm_ids)] <- probe_ids[is.na(norm_ids)]
  
  if (!anyDuplicated(norm_ids)) {
    rownames(mat) <- norm_ids
    return(mat)
  }
  
  if (fun == "mean") {
    collapsed <- rowsum(mat, group = norm_ids, reorder = FALSE) /
      as.vector(table(norm_ids)[rownames(rowsum(mat, group = norm_ids, reorder = FALSE))])
    return(as.matrix(collapsed))
  } else {
    split_idx <- split(seq_along(norm_ids), norm_ids)
    out <- do.call(rbind, lapply(split_idx, function(idx) {
      if (length(idx) == 1) {
        mat[idx, , drop = FALSE]
      } else {
        apply(mat[idx, , drop = FALSE], 2, median, na.rm = TRUE) |>
          matrix(nrow = 1, dimnames = list(names(split_idx)[which(names(split_idx) == names(split_idx)[1])], colnames(mat)))
      }
    }))
    rownames(out) <- names(split_idx)
    return(as.matrix(out))
  }
}

save_outputs <- function(beta_mat, mval_mat, sample_meta, out_dir) {
  # Beta matrix with CpG IDs preserved
  beta_df_out <- as.data.frame(beta_mat, check.names = FALSE)
  beta_df_out <- cbind(CpG = rownames(beta_df_out), beta_df_out)
  
  saveRDS(beta_df_out, file.path(out_dir, "salas_beta_matrix.rds"))
  fwrite(beta_df_out, file.path(out_dir, "salas_beta_matrix.csv"))
  
  # M-value matrix with CpG IDs preserved
  mval_df_out <- as.data.frame(mval_mat, check.names = FALSE)
  mval_df_out <- cbind(CpG = rownames(mval_df_out), mval_df_out)
  
  saveRDS(mval_df_out, file.path(out_dir, "salas_mval_matrix.rds"))
  fwrite(mval_df_out, file.path(out_dir, "salas_mval_matrix.csv"))
  
  # Sample metadata
  saveRDS(sample_meta, file.path(out_dir, "salas_sample_metadata.rds"))
  fwrite(sample_meta, file.path(out_dir, "salas_sample_metadata.csv"))
}

# ------------------------------------------------
# Main
# ------------------------------------------------

message2("Reading metadata: %s", META_CSV)
meta <- read_metadata(META_CSV)
message2("Metadata rows: %d", nrow(meta))

message2("Finding IDAT pairs under: %s", IDAT_ROOT)
idat_df <- find_idat_pairs(IDAT_ROOT)
message2("Found %d complete IDAT pairs", nrow(idat_df))

message2("Matching IDAT pairs to metadata")
matched_df <- match_idats_to_metadata(idat_df, meta)

fwrite(matched_df, file.path(OUT_DIR, "idat_metadata_match_table.csv"))

message2("Matched %d / %d IDAT pairs", sum(matched_df$matched), nrow(matched_df))

# Keep only matched samples with non-missing cell type
use_df <- matched_df |>
  filter(matched, !is.na(cell_type), cell_type != "")

if (nrow(use_df) == 0) {
  stop("No matched IDATs with cell_type were found.")
}

message2("Using %d samples with assigned cell_type", nrow(use_df))
print(table(use_df$cell_type))

# Optional: keep only your target cell types
target_cell_types <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                       "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")

use_df <- use_df |>
  filter(cell_type %in% target_cell_types)

message2("Retained %d samples in target 12-cell schema", nrow(use_df))
print(table(use_df$cell_type))

# ------------------------------------------------
# Process IDATs with SeSAMe
# ------------------------------------------------
beta_list <- list()
mval_list <- list()
sample_rows <- vector("list", nrow(use_df))

for (i in seq_len(nrow(use_df))) {
  row <- use_df[i, ]
  base <- row$basename
  
  sample_name <- if (!is.na(row$gsm) && nzchar(row$gsm)) {
    row$gsm
  } else {
    row$sample_stub
  }
  
  message2("[%d/%d] Processing %s | %s | %s", i, nrow(use_df), sample_name, row$gse, row$cell_type)
  
  res <- prep_sesame_one(base, prep = "QCDPB")
  
  beta_list[[sample_name]] <- res$beta
  mval_list[[sample_name]] <- res$mval
  
  sample_rows[[i]] <- data.frame(
    sample_name = sample_name,
    gsm = row$gsm,
    gse = row$gse,
    cell_type = row$cell_type,
    basename = row$basename,
    sample_stub = row$sample_stub,
    sample_title = row$Sample_title %||% NA_character_,
    source_name = row$Sample_source_name_ch1 %||% NA_character_,
    characteristics_text = row$characteristics_text %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

# helper for null-coalescing
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

sample_meta <- bind_rows(sample_rows)

# ------------------------------------------------
# Align probes across all samples
# ------------------------------------------------
message2("Aligning probes across samples")

common_beta_ids <- Reduce(intersect, lapply(beta_list, names))
common_mval_ids <- Reduce(intersect, lapply(mval_list, names))
common_ids <- intersect(common_beta_ids, common_mval_ids)

if (length(common_ids) == 0) {
  stop("No common probes found across all samples.")
}

message2("Common probes across all samples: %d", length(common_ids))

beta_mat <- do.call(cbind, lapply(beta_list, function(x) x[common_ids]))
mval_mat <- do.call(cbind, lapply(mval_list, function(x) x[common_ids]))

colnames(beta_mat) <- names(beta_list)
colnames(mval_mat) <- names(mval_list)
rownames(beta_mat) <- common_ids
rownames(mval_mat) <- common_ids

# ------------------------------------------------
# Optional: collapse EPICv2 replicate suffixes
# ------------------------------------------------
message2("Collapsing duplicate/replicate probe suffixes by mean")
beta_mat <- collapse_duplicate_probes(beta_mat, fun = "mean")
mval_mat <- collapse_duplicate_probes(mval_mat, fun = "mean")

message2("Final beta matrix dimensions: %d x %d", nrow(beta_mat), ncol(beta_mat))

# ------------------------------------------------
# Basic QC summaries
# ------------------------------------------------
na_per_sample <- colSums(is.na(beta_mat))
na_frac_per_sample <- na_per_sample / nrow(beta_mat)

qc_df <- data.frame(
  sample_name = colnames(beta_mat),
  n_na = na_per_sample,
  frac_na = na_frac_per_sample,
  stringsAsFactors = FALSE
) |>
  left_join(sample_meta, by = "sample_name")

fwrite(qc_df, file.path(OUT_DIR, "sample_qc_summary.csv"))

message2("NA fraction summary:")
print(summary(na_frac_per_sample))

message2("Worst samples by NA fraction:")
print(head(qc_df[order(qc_df$frac_na, decreasing = TRUE), c("sample_name", "cell_type", "frac_na")], 20))

# ------------------------------------------------
# Save outputs
# ------------------------------------------------
save_outputs(beta_mat, mval_mat, sample_meta, OUT_DIR)

# ------------------------------------------------
# Build per-cell-type centroids (useful for deconvolution library)
# ------------------------------------------------
message2("Computing per-cell-type centroids")

centroid_mean <- sapply(target_cell_types, function(ct) {
  cols <- sample_meta$sample_name[sample_meta$cell_type == ct]
  if (length(cols) == 0) {
    rep(NA_real_, nrow(beta_mat))
  } else if (length(cols) == 1) {
    beta_mat[, cols]
  } else {
    rowMeans(beta_mat[, cols, drop = FALSE], na.rm = TRUE)
  }
})

centroid_median <- sapply(target_cell_types, function(ct) {
  cols <- sample_meta$sample_name[sample_meta$cell_type == ct]
  if (length(cols) == 0) {
    rep(NA_real_, nrow(beta_mat))
  } else if (length(cols) == 1) {
    beta_mat[, cols]
  } else {
    apply(beta_mat[, cols, drop = FALSE], 1, median, na.rm = TRUE)
  }
})

rownames(centroid_mean) <- rownames(beta_mat)
rownames(centroid_median) <- rownames(beta_mat)

saveRDS(centroid_mean, file.path(OUT_DIR, "salas_centroid_mean_beta.rds"))
saveRDS(centroid_median, file.path(OUT_DIR, "salas_centroid_median_beta.rds"))

fwrite(as.data.frame(centroid_mean, check.names = FALSE),
       file.path(OUT_DIR, "salas_centroid_mean_beta.csv"))
fwrite(as.data.frame(centroid_median, check.names = FALSE),
       file.path(OUT_DIR, "salas_centroid_median_beta.csv"))

# ------------------------------------------------
# Save a Python-friendly metadata table
# ------------------------------------------------
python_meta <- sample_meta |>
  select(sample_name, gsm, gse, cell_type)

fwrite(python_meta, file.path(OUT_DIR, "salas_python_sample_to_celltype.csv"))

message2("Done.")
message2("Outputs written to: %s", OUT_DIR)

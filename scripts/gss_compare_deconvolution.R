#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sesame)
  library(EpiDISH)
  library(data.table)
})

# ------------------------------------------------
# Inputs
# ------------------------------------------------
truth_file <- "gss/salas_mix_proportions.csv"

normalize_mix_name <- function(x) {
  # "mem mix15" -> "memmix15"
  gsub("\\s+", "", trimws(x))
}

extract_sentrix_id <- function(x) {
  # From:
  # /path/GSM5121368_204361720136_R01C01
  # to:
  # 204361720136_R01C01
  
  b <- basename(x)
  sub("^GSM[0-9]+_", "", b)
}

IDAT_DIRS <- c(
  "gss/idat/GSE182379",
  "gss/idat/GSE167998"
)

assaf_ref_file <- "gss/assaf_reference_matrix.csv"
gss_ref_file   <- "gss/gss_reference_matrix.csv"

out_dir <- "results/deconvolution/comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helpers
# ------------------------------------------------
extract_gsm <- function(x) {
  sub("^([^_]+)_.*$", "\\1", basename(x))
}

load_reference_csv <- function(path) {
  ref_df <- read.csv(path, check.names = FALSE)
  rownames(ref_df) <- ref_df[[1]]
  ref_df[[1]] <- NULL
  
  ref <- as.matrix(ref_df)
  mode(ref) <- "numeric"
  
  # Ensure cg-only rownames
  rownames(ref) <- sub("_.*$", "", rownames(ref))
  ref
}

collapse_beta_to_cg <- function(beta_mat) {
  cg_ids <- sub("_.*$", "", rownames(beta_mat))
  
  # If no duplicates, just rename rows
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

calc_metrics <- function(est, truth) {
  stopifnot(identical(rownames(est), rownames(truth)))
  stopifnot(identical(colnames(est), colnames(truth)))
  
  pcc <- sapply(colnames(est), function(ct) {
    suppressWarnings(cor(est[, ct], truth[, ct], use = "pairwise.complete.obs"))
  })
  
  rmse <- sapply(colnames(est), function(ct) {
    sqrt(mean((est[, ct] - truth[, ct])^2, na.rm = TRUE))
  })
  
  global_rmse <- sqrt(mean((as.matrix(est) - as.matrix(truth))^2, na.rm = TRUE))
  
  list(
    per_celltype = data.frame(
      CellType = colnames(est),
      PCC = pcc,
      RMSE = rmse,
      row.names = NULL
    ),
    global_rmse = global_rmse
  )
}

run_deconv_with_reference <- function(beta_mix, ref_cg, max_na_frac = 0.20) {
  beta_cg_mat <- collapse_beta_to_cg(beta_mix)
  
  common_cg <- intersect(rownames(ref_cg), rownames(beta_cg_mat))
  
  if (length(common_cg) < 50) {
    stop("Too few shared probes: ", length(common_cg))
  }
  
  beta_sub <- beta_cg_mat[common_cg, , drop = FALSE]
  ref_sub  <- ref_cg[common_cg, , drop = FALSE]
  
  stopifnot(identical(rownames(beta_sub), rownames(ref_sub)))
  
  na_frac <- colSums(is.na(beta_sub)) / nrow(beta_sub)
  keep_samples <- na_frac <= max_na_frac
  
  if (!all(keep_samples)) {
    message("Dropping samples with high NA fraction:")
    print(names(keep_samples)[!keep_samples])
  }
  
  beta_sub <- beta_sub[, keep_samples, drop = FALSE]
  
  # Remove probes with remaining NA after sample filtering
  keep_probes <- complete.cases(beta_sub) & complete.cases(ref_sub)
  beta_sub <- beta_sub[keep_probes, , drop = FALSE]
  ref_sub  <- ref_sub[keep_probes, , drop = FALSE]
  
  message("Shared probes used: ", nrow(beta_sub))
  message("Samples used: ", ncol(beta_sub))
  
  deconv <- epidish(
    beta.m = beta_sub,
    ref.m  = ref_sub,
    method = "RPC"
  )
  
  list(
    est = deconv$estF,
    beta_sub = beta_sub,
    ref_sub = ref_sub,
    kept_samples = colnames(beta_sub)
  )
}

# ------------------------------------------------
# Discover IDATs
# ------------------------------------------------
message("Building beta matrix from IDAT files...")

idat_files <- unlist(lapply(IDAT_DIRS, function(d) {
  list.files(d, pattern = "_Red.idat.gz$", full.names = TRUE)
}))

basenames_full <- sub("_Red.idat.gz$", "", idat_files)

message("Found ", length(basenames_full), " IDAT pairs across all directories")

idat_df <- data.frame(
  basename = basenames_full,
  gsm = extract_gsm(basenames_full),
  sentrix_id = extract_sentrix_id(basenames_full),
  stringsAsFactors = FALSE
)

truth_dt <- fread(truth_file)

required_cols <- c("Mix Name", "Sentrix ID")
missing_cols <- setdiff(required_cols, colnames(truth_dt))
if (length(missing_cols) > 0) {
  stop("salas_mix_proportions.csv is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

truth_dt[, Sample_title := normalize_mix_name(`Mix Name`)]
truth_dt[, sentrix_id := trimws(`Sentrix ID`)]

meta_merged <- merge(
  idat_df,
  truth_dt[, .(sentrix_id, Sample_title)],
  by = "sentrix_id"
)

message("Matched ", nrow(meta_merged), " samples to IDATs using Sentrix ID")

if (nrow(meta_merged) == 0) {
  stop("No IDAT files matched salas_mix_proportions.csv by Sentrix ID")
}

message("Matched samples:")
print(meta_merged[, c("gsm", "sentrix_id", "Sample_title")])

# ------------------------------------------------
# Read + preprocess IDATs
# ------------------------------------------------
beta_list <- list()

for (i in seq_len(nrow(meta_merged))) {
  if (i %% 10 == 0 || i == 1) {
    message(sprintf("Processing %d / %d", i, nrow(meta_merged)))
  }
  
  base <- meta_merged$basename[i]
  sample_name <- meta_merged$Sample_title[i]
  
  s <- readIDATpair(base)
  s <- prepSesame(s, "QCDPB")
  beta <- getBetas(s)
  
  beta_list[[sample_name]] <- beta
}

common_cpgs <- Reduce(intersect, lapply(beta_list, names))

beta_mat <- do.call(cbind, lapply(beta_list, function(b) b[common_cpgs]))
colnames(beta_mat) <- names(beta_list)
rownames(beta_mat) <- common_cpgs

message("Beta matrix dims: ", paste(dim(beta_mat), collapse = " x "))
message("Beta sample names:")
print(colnames(beta_mat))

# ------------------------------------------------
# Load mixture truth table
# ------------------------------------------------
message("Loading corrected Salas mixture truth table...")

truth_dt <- fread(truth_file)

truth_dt[, Sample_title := normalize_mix_name(`Mix Name`)]

celltype_cols <- setdiff(
  colnames(truth_dt),
  c("Mix Name", "Sentrix ID", "Sample_title")
)

truth_df <- as.data.frame(
  truth_dt[, c("Sample_title", celltype_cols), with = FALSE],
  check.names = FALSE
)

rownames(truth_df) <- truth_df$Sample_title
truth_df$Sample_title <- NULL

truth_df[] <- lapply(truth_df, as.numeric)

if (max(as.matrix(truth_df), na.rm = TRUE) > 1.5) {
  message("Converting truth table from percent to fraction")
  truth_df <- truth_df / 100
}

message("Truth table dims: ", paste(dim(truth_df), collapse = " x "))
message("Truth sample names:")
print(rownames(truth_df))

message("Truth range:")
print(range(as.matrix(truth_df), na.rm = TRUE))

message("Truth row sums:")
print(rowSums(truth_df))

# ------------------------------------------------
# Align beta and truth
# ------------------------------------------------
common_samples <- intersect(colnames(beta_mat), rownames(truth_df))

message("Common samples: ", length(common_samples))
print(common_samples)

if (length(common_samples) == 0) {
  stop("No overlap between beta matrix samples and truth table")
}

missing_truth <- setdiff(colnames(beta_mat), rownames(truth_df))
missing_beta  <- setdiff(rownames(truth_df), colnames(beta_mat))

if (length(missing_truth) > 0) {
  message("Samples in beta matrix but not in truth table:")
  print(missing_truth)
}

if (length(missing_beta) > 0) {
  message("Samples in truth table but not in beta matrix:")
  print(missing_beta)
}

beta_mix <- beta_mat[, common_samples, drop = FALSE]
truth_df <- truth_df[common_samples, , drop = FALSE]

stopifnot(identical(colnames(beta_mix), rownames(truth_df)))

message("Aligned samples: ", length(common_samples))

# Remove CpGs with NA before deconvolution
message("Before NA filtering: ", paste(dim(beta_mix), collapse = " x "))
beta_mix <- beta_mix[complete.cases(beta_mix), , drop = FALSE]
message("After NA filtering: ", paste(dim(beta_mix), collapse = " x "))

# ------------------------------------------------
# Load references
# ------------------------------------------------
assaf_ref <- load_reference_csv(assaf_ref_file)
gss_ref   <- load_reference_csv(gss_ref_file)

cat("Assaf ref dims:", dim(assaf_ref), "\n")
cat("GSS ref dims:", dim(gss_ref), "\n")

message("Assaf ref range:")
print(range(assaf_ref, na.rm = TRUE))
message("GSS ref range:")
print(range(gss_ref, na.rm = TRUE))
message("Beta range:")
print(range(beta_mix, na.rm = TRUE))

# ------------------------------------------------
# Harmonize cell-type columns
# ------------------------------------------------
common_celltypes <- Reduce(intersect, list(
  colnames(assaf_ref),
  colnames(gss_ref),
  colnames(truth_df)
))

if (length(common_celltypes) < 2) {
  stop("Too few shared cell types across references and truth.")
}

assaf_ref <- assaf_ref[, common_celltypes, drop = FALSE]
gss_ref   <- gss_ref[, common_celltypes, drop = FALSE]
truth_df  <- truth_df[, common_celltypes, drop = FALSE]

cat("Shared cell types:", paste(common_celltypes, collapse = ", "), "\n")

# ------------------------------------------------
# Run deconvolution
# ------------------------------------------------
assaf_run <- run_deconv_with_reference(beta_mix, assaf_ref)
gss_run   <- run_deconv_with_reference(beta_mix, gss_ref)

# ------------------------------------------------
# Align estimates and truth
# ------------------------------------------------
est_assaf <- assaf_run$est[, common_celltypes, drop = FALSE]
est_gss   <- gss_run$est[, common_celltypes, drop = FALSE]

truth_assaf <- truth_df[rownames(est_assaf), common_celltypes, drop = FALSE]
truth_gss   <- truth_df[rownames(est_gss), common_celltypes, drop = FALSE]

stopifnot(identical(rownames(est_assaf), rownames(truth_assaf)))
stopifnot(identical(rownames(est_gss), rownames(truth_gss)))

# ------------------------------------------------
# Metrics
# ------------------------------------------------
assaf_metrics <- calc_metrics(est_assaf, truth_assaf)
gss_metrics   <- calc_metrics(est_gss, truth_gss)

metrics_df <- merge(
  assaf_metrics$per_celltype,
  gss_metrics$per_celltype,
  by = "CellType",
  suffixes = c("_Assaf", "_GSS")
)

metrics_df$Delta_RMSE <- metrics_df$RMSE_GSS - metrics_df$RMSE_Assaf
metrics_df$Delta_PCC  <- metrics_df$PCC_GSS - metrics_df$PCC_Assaf

global_df <- data.frame(
  Reference = c("Assaf", "GSS"),
  Global_RMSE = c(assaf_metrics$global_rmse, gss_metrics$global_rmse),
  Samples_Used = c(nrow(est_assaf), nrow(est_gss)),
  Shared_Probes = c(nrow(assaf_run$ref_sub), nrow(gss_run$ref_sub))
)

# ------------------------------------------------
# Save outputs
# ------------------------------------------------
write.csv(metrics_df,
          file.path(out_dir, "per_celltype_metrics_comparison_new.csv"),
          row.names = FALSE)

write.csv(global_df,
          file.path(out_dir, "global_metrics_comparison_new.csv"),
          row.names = FALSE)

write.csv(cbind(Sample = rownames(est_assaf), est_assaf),
          file.path(out_dir, "estimated_fractions_assaf_new.csv"),
          row.names = FALSE)

write.csv(cbind(Sample = rownames(est_gss), est_gss),
          file.path(out_dir, "estimated_fractions_gss_new.csv"),
          row.names = FALSE)

write.csv(cbind(Sample = rownames(truth_df), truth_df),
          file.path(out_dir, "truth_fractions_aligned_new.csv"),
          row.names = FALSE)

cat("\nGlobal metrics:\n")
print(global_df)

cat("\nPer-cell-type metrics:\n")
print(metrics_df)

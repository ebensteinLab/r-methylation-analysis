suppressPackageStartupMessages({
  library(matrixStats)
  library(dplyr)
  library(sva)
})

compute_marker_fractions <- function(beta_mat, probe_map) {
  # beta_mat: probes x samples
  # Returns: samples x cell_types
  
  cell_types <- names(probe_map)
  
  # Activity matrix: 1 - beta
  act <- 1 - beta_mat
  
  scores <- sapply(cell_types, function(ct) {
    probes <- probe_map[[ct]]
    # median across probes for each sample
    colMedians(act[probes, , drop = FALSE], na.rm = TRUE)
  })
  scores[is.nan(scores)] <- NA
  
  print(head(scores))
  print(apply(scores, 2, summary))
  print(apply(scores, 2, sd, na.rm=TRUE))
  
  # scores: samples x cell_types
  # convert to pseudo-fractions per sample
  fracs <- t(apply(scores, 1, function(x) {
    x[x < 0] <- 0
    x / sum(x)
  }))
  colnames(fracs) <- cell_types
  fracs
}

probe_map <- list(
  B    = c("cg02212339_BC11","cg19045191_BC21","cg17232476_BC21"),
  CD8  = c("cg00219921_TC21"),
  Eos  = c("cg02803925_TC11","cg05736642_TC21"),
  Mono = c("cg24788483_TC11"),
  Neut = c("cg13595556_TC21","cg09859659_TC21","cg25600606_TC21"),
  T    = c("cg14001486_BC21","cg09088496_TC21"),
  Treg = c("cg02033323_BC21")
)

# Inputs
gen_beta <- readRDS("results/processed/beta_matrix_sesame.rds")
targets  <- readRDS("results/processed/targets_with_sesame.rds")
batch <- targets$Sentrix_Id
gen_beta <- ComBat(dat = gen_beta, batch = batch)

cf_beta  <- readRDS("results/processed/cf_beta_matrix_sesame.rds")
cf_targets  <- readRDS("results/processed/cf_targets_with_sesame.rds")
cf_batch <- cf_targets$Sentrix_Id
cf_beta  <- ComBat(dat = cf_beta,  batch = cf_batch)

all_probes <- unique(unlist(probe_map))
missing_gen <- setdiff(all_probes, rownames(gen_beta))
missing_cf  <- setdiff(all_probes, rownames(cf_beta))

cat("Missing in gen_beta:", paste(missing_gen, collapse=", "), "\n")
cat("Missing in cf_beta :", paste(missing_cf, collapse=", "), "\n")
stopifnot(length(missing_gen) == 0, length(missing_cf) == 0)

gen_marker_fracs <- compute_marker_fractions(gen_beta[all_probes, , drop=FALSE], probe_map)
cf_marker_fracs  <- compute_marker_fractions(cf_beta[all_probes, , drop=FALSE], probe_map)

# Make sure Patient exists
stopifnot("Patient" %in% colnames(targets))

# Trim whitespace just in case
colnames(gen_beta) <- trimws(colnames(gen_beta))
targets$Patient <- trimws(targets$Patient)

idx <- match(colnames(gen_beta), targets$Patient)

# Debug: how many didn't match?
cat("Unmatched:", sum(is.na(idx)), "of", length(idx), "\n")
if (anyNA(idx)) {
  print(head(colnames(gen_beta)[is.na(idx)], 20))
  print(head(targets$Patient, 20))
}

# Reorder
targets2 <- targets[idx, ]

# Now this should be TRUE
stopifnot(identical(colnames(gen_beta), targets2$Patient))

colnames(cf_beta) <- trimws(colnames(cf_beta))
idx_cf <- match(colnames(cf_beta), targets$Patient)
cf_names <- trimws(colnames(cf_beta))
pat_names <- trimws(targets$Patient)

idx_cf <- match(cf_names, pat_names)

cat("Unmatched cf:", sum(is.na(idx_cf)), "of", length(idx_cf), "\n")
targets_cf <- targets[idx_cf, ]
stopifnot(identical(colnames(cf_beta), targets_cf$Patient))

gen_out <- cbind(targets2, as.data.frame(gen_marker_fracs))
cf_out <- cbind(targets_cf, as.data.frame(cf_marker_fracs))

write.csv(gen_out, "results/deconvolution/marker_pseudofractions_genomic.csv", row.names = FALSE)
write.csv(cf_out,  "results/deconvolution/marker_pseudofractions_cf.csv", row.names = FALSE)

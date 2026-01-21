#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

library(ExperimentHub)
library(minfi)
library(limma)
library(IDOL)
library(sesameData)
library(rtracklayer)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

sesameDataCache()

# ------------------------------------------------------------
# Load FlowSorted Blood EPIC reference
# ------------------------------------------------------------
eh <- ExperimentHub()
ref_rg <- eh[["EH1136"]]   # EPIC blood reference

ref_mset <- preprocessNoob(ref_rg)
ref_beta <- getBeta(ref_mset)
ref_mval <- getM(ref_mset)
celltypes <- factor(pData(ref_mset)$CellType)

rm(ref_rg, ref_mset)
gc()

message("Reference samples: ", length(celltypes))
message("Cell types: ", paste(levels(celltypes), collapse = ", "))

ref_beta_centroids <- sapply(levels(celltypes), function(ct) {
  rowMeans(ref_beta[, celltypes == ct, drop = FALSE], na.rm = TRUE)
})
ref_beta_centroids <- as.matrix(ref_beta_centroids)

saveRDS(ref_beta_centroids, "results/processed/ref_beta_epic_centroids.rds")
rm(ref_beta_centroids)
gc()

# ------------------------------------------------------------
# Load SeSAMe bulk M-values
# ------------------------------------------------------------
mval_bulk <- readRDS("results/processed/mval_matrix_sesame_batch_corrected.rds")
targets <- readRDS("results/processed/targets_merged.rds")

# ------------------------------------------------------------
# Harmonize SeSAMe and EPIC probe IDs
# ------------------------------------------------------------

# Map EPICv2 probe IDs -> EPIC v1 probe IDs
addr <- sesameDataGet("EPICv2.address")
epicv2_gr <- addr$hg38
names(epicv2_gr) <- names(addr$hg38)

epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

epic_gr <- GRanges(
  seqnames = epic_anno$chr,
  ranges   = IRanges(epic_anno$pos, width=1),
  strand   = "*",
  probe    = rownames(epic_anno)
)

rm(epic_anno, addr)

chain <- import.chain("hg19ToHg38.over.chain")

epic_gr_hg38 <- liftOver(epic_gr, chain)
epic_gr_hg38 <- unlist(epic_gr_hg38)

# --- EPIC (hg19) -> hg38 via liftOver, preserving EPIC probe IDs ---
epic_gr_hg38_list <- liftOver(epic_gr, chain)

# Keep only probes that map to exactly one hg38 locus (avoid multi-mappers)
one_map <- lengths(epic_gr_hg38_list) == 1
epic_gr_hg38 <- unlist(epic_gr_hg38_list[one_map], use.names = FALSE)

# Restore EPIC probe IDs as names after unlist
names(epic_gr_hg38) <- epic_gr$probe[one_map]

# Ensure widths match EPICv2 (EPICv2 ranges are width=2 in your object)
epic_gr_hg38 <- resize(epic_gr_hg38, width = 2, fix = "start")

# --- IMPORTANT: overlap in the SAME genome build (hg38 vs hg38) ---
hits <- findOverlaps(epic_gr_hg38, epicv2_gr, maxgap = 0)

# EPICv2 probe IDs
epicv2_ids <- names(epicv2_gr)[subjectHits(hits)]
# EPIC probe IDs (from names we set)
epic_ids_mapped <- names(epic_gr_hg38)[queryHits(hits)]

rm(hits, epic_gr_hg38, epic_gr_hg38_list, one_map, chain)

# Build lookup: EPICv2 -> EPIC
# If multiple EPIC map to same EPICv2, keep the first (you can change policy)
epicv2_to_epic <- epic_ids_mapped
names(epicv2_to_epic) <- epicv2_ids
epicv2_to_epic <- epicv2_to_epic[!duplicated(names(epicv2_to_epic))]
saveRDS(epicv2_to_epic, "results/processed/epicv2_to_epic_map.rds")

rm(epicv2_gr, epic_ids_mapped)

message("Mapped EPICv2 IDs: ", length(epicv2_to_epic))

mapped_fraction <- mean(rownames(mval_bulk) %in% names(epicv2_to_epic))
message("Fraction of bulk probes that have a mapping: ", round(mapped_fraction, 3))

# Convert SeSAMe rownames
old_ids <- rownames(mval_bulk)
new_ids <- epicv2_to_epic[old_ids]

rownames(mval_bulk) <- new_ids
mval_bulk <- mval_bulk[!is.na(rownames(mval_bulk)), , drop = FALSE]

epic_gr <- sesameData_getManifestGRanges("EPIC")
epic_ids <- names(epic_gr)   # EPIC cg IDs
message("EPIC IDs: ", length(epic_ids))
message("Mapped EPICv2 IDs: ", length(epicv2_to_epic))
message("N bulk probes: ", nrow(mval_bulk))
message("N mapped: ", sum(!is.na(new_ids)))

common <- intersect(rownames(ref_mval), rownames(mval_bulk))
common <- intersect(common, epic_ids)

message("Shared CpGs: ", length(common))
stopifnot(length(common) > 300000)

rm(epic_gr, epic_ids, epicv2_to_epic, old_ids, new_ids)

ref_beta <- ref_beta[common, ]

saveRDS(ref_beta, "results/processed/ref_beta_epic.rds")

rm(common)

gc()


# ------------------------------------------------------------
# Build candidate DMR finder (IDOL internal)
# ------------------------------------------------------------
message("Creating candFinder")
covars <- data.frame(CellType = celltypes, dummy = seq_along(celltypes))
rownames(covars) <- colnames(ref_beta)
candFinder <- CandidateDMRFinder.v2(
  cellTypes       = celltypes,
  referenceBetas  = ref_beta,
  referenceCovars = covars,
  M = 150,
  equal.variance = FALSE
)

# ------------------------------------------------------------
# Patch candFinder to make it compatible with IDOLoptimize
# ------------------------------------------------------------

coef <- candFinder$coefEsts
labs <- colnames(coef)

# unique cell types
u <- unique(labs)

# collapse coefEsts by cell type (mean across samples)
coef_collapsed <- sapply(u, function(ct) {
  rowMeans(coef[, labs == ct, drop = FALSE])
})

# Replace coefEsts with collapsed version
candFinder$coefEsts <- coef_collapsed
message(" ", colnames(candFinder$coefEsts))

# ------------------------------------------------------------
# IDOL optimization
# ------------------------------------------------------------

idol_classes <- colnames(candFinder$coefEsts)

onehot <- model.matrix(~ 0 + celltypes)
colnames(onehot) <- levels(celltypes)

trainingCovariates <- onehot[, idol_classes, drop = FALSE]
rownames(trainingCovariates) <- colnames(ref_beta)

stopifnot(all(colnames(trainingCovariates) == idol_classes))

rm(onehot, idol_classes, celltypes)
gc()

message("Running IDOLoptimize")
idol_res <- IDOLoptimize(
  candDMRFinderObject = candFinder,
  trainingBetas       = ref_beta,
  trainingCovariates  = trainingCovariates,
  libSize             = 300,
  maxIt               = 500,
  numCores            = 4
)

rm(candFinder, ref_beta, trainingCovariates)
gc()

idol_probes <- idol_res[[ "IDOL Optimized Library" ]]

message("IDOL selected CpGs: ", length(idol_probes))

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
saveRDS(idol_probes, "results/processed/idol_cpgs.rds")
saveRDS(idol_res, "results/processed/idol_model.rds")

rm(idol_probes, idol_res)
gc()

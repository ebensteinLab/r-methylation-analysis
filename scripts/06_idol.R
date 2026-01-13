# -----------------------------------------------------------
# Install IDOL if needed
# -----------------------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("immunomethylomics/IDOL")

BiocManager::install(c(
  "EpiDISH",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "limma",
  "GenomicRanges",
  "DMRcate"
))

library(ExperimentHub)
library(minfi)
library(limma)
library(IDOL)
library(sesame)
library(sesameData)
library(DMRcate)
library(rtracklayer)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

set.seed(1)

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

message("Reference samples: ", length(celltypes))
message("Cell types: ", paste(levels(celltypes), collapse = ", "))

# ------------------------------------------------------------
# Load SeSAMe bulk M-values
# ------------------------------------------------------------
mval_bulk <- readRDS("results/processed/mval_matrix_sesame_batch_corrected.rds")

# ------------------------------------------------------------
# Harmonize SeSAMe and EPIC probe IDs
# ------------------------------------------------------------
sesameDataCache()

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

# Build lookup: EPICv2 -> EPIC
# If multiple EPIC map to same EPICv2, keep the first (you can change policy)
epicv2_to_epic <- epic_ids_mapped
names(epicv2_to_epic) <- epicv2_ids
epicv2_to_epic <- epicv2_to_epic[!duplicated(names(epicv2_to_epic))]

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

ref_mval <- ref_mval[common, ]
ref_beta <- ref_beta[common, ]
mval_bulk <- mval_bulk[common, ]

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
# Run IDOL optimization
# ------------------------------------------------------------

idol_classes <- colnames(candFinder$coefEsts)

onehot <- model.matrix(~ 0 + celltypes)
colnames(onehot) <- levels(celltypes)

trainingCovariates <- onehot[, idol_classes, drop = FALSE]
rownames(trainingCovariates) <- colnames(ref_beta)

stopifnot(all(colnames(trainingCovariates) == idol_classes))

message("Running IDOLoptimize")
idol_res <- IDOLoptimize(
  candDMRFinderObject = candFinder,
  trainingBetas       = ref_beta,
  trainingCovariates  = trainingCovariates,
  libSize             = 300,
  maxIt               = 500,
  numCores            = 4
)

idol_probes <- idol_res[[ "IDOL Optimized Library" ]]

message("IDOL selected CpGs: ", length(idol_probes))

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
saveRDS(idol_probes, "results/processed/idol_cpgs.rds")
saveRDS(idol_res,    "results/processed/idol_model.rds")

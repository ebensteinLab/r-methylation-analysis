# install.packages("Matrix")
# 
# library(Matrix)
# 
# install.packages("BiocManager")
# BiocManager::install("genefilter", force = TRUE)
# 
# library(genefilter)

beta_df <- read.csv("gss/processed/salas_beta_matrix.csv", check.names = FALSE)

rownames(beta_df) <- beta_df$CpG
beta_df$CpG <- NULL

beta_mat <- as.matrix(beta_df)
mode(beta_mat) <- "numeric"

beta_mat_noNA <- beta_mat[complete.cases(beta_mat), ]

meta <- read.csv("gss/processed/salas_sample_metadata.csv")

all(colnames(beta_mat_noNA) == meta$sample_name)

cell.types.info <- data.frame(
  CellType = meta$cell_type
)

panel_cell_types <- unique(cell.types.info$CellType)

source("~/R/projects/r-methylation-analysis/scripts/ConstructDNAmPanel.R")

res <- ConstructDNAmPanel(
  dnam.matrix = beta_mat_noNA,
  cell.types.info = cell.types.info,
  panel.cell.types = panel_cell_types,
  n.probes = 50,
  fdr.threshold = 0.05,
  summary.method = "mean",
  equal.variance = FALSE
)

saveRDS(res, "gss/salas_panel_full.rds")

ref <- res$ReferencePanel

saveRDS(ref, "gss/salas_reference_matrix.rds")

# also export to CSV for Python
ref_df <- cbind(CpG = rownames(ref), as.data.frame(ref, check.names = FALSE))
data.table::fwrite(ref_df, "gss/salas_reference_matrix.csv")

saveRDS(res$DetailedResults, "gss/salas_marker_details.rds")


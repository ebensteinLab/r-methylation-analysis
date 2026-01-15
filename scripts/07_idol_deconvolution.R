#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(EpiDISH)
})

message("Loading inputs...")

# ------------------------------------------------
# Load data
# ------------------------------------------------
beta <- readRDS("results/processed/beta_matrix_sesame_batch_corrected.rds")
targets <- readRDS("results/processed/targets_with_sesame.rds")
epicv2_to_epic <- readRDS("results/processed/epicv2_to_epic_map.rds")
idol_probes <- readRDS("results/processed/idol_cpgs_plasma.rds")

old_ids <- rownames(beta)
new_ids <- epicv2_to_epic[old_ids]

# Drop probes without mapping
keep <- !is.na(new_ids)

beta <- beta[keep, , drop = FALSE]
rownames(beta) <- new_ids[keep]

# ------------------------------------------------
# Identify plasma samples
# ------------------------------------------------
is_plasma <- grepl("(MG|SMM|MM)$", targets$Patient)

stopifnot(
  ncol(beta) == nrow(targets),
  length(is_plasma) == ncol(beta)
)

message("Plasma samples: ", sum(is_plasma))

# ------------------------------------------------
# Subset bulk betas
# ------------------------------------------------
beta_plasma <- beta[, is_plasma, drop = FALSE]

# ------------------------------------------------
# Restrict to IDOL CpGs
# ------------------------------------------------
common_idol <- Reduce(
  intersect,
  list(
    idol_probes,
    rownames(beta_plasma),
    rownames(ref_beta_centroids)
  )
)

message("IDOL CpGs used: ", length(common_idol))
stopifnot(length(common_idol) > 200)

# ------------------------------------------------
# Run EpiDISH (RPC)
# ------------------------------------------------
message("Running EpiDISH deconvolution (RPC)...")

# bulk: CpGs x bulk samples
beta_plasma_idol <- beta_plasma[common_idol, , drop = FALSE]

# reference: CpGs x celltypes  (centroids)
ref_beta_centroids <- readRDS("results/processed/ref_beta_epic_centroids.rds")
ref_beta_idol <- ref_beta_centroids[common_idol, , drop = FALSE]

deconv <- epidish(beta.m = beta_plasma_idol, ref.m = ref_beta_idol, method = "RPC")

fractions <- deconv$estF

# ------------------------------------------------
# Save output
# ------------------------------------------------
out_file <- "results/processed/plasma_cell_fractions_idol.rds"
saveRDS(fractions, out_file)

message("Saved plasma deconvolution results to:")
message(out_file)

plasma_targets <- targets[is_plasma, ]

# align order with fractions (VERY IMPORTANT)
plasma_targets <- plasma_targets[
  match(rownames(fractions), plasma_targets$sample_name),
]

# derive disease
plasma_targets$Disease <- substr(plasma_targets$Patient, 6, nchar(plasma_targets$Patient))

df <- cbind(plasma_targets[, c("sample_name", "Patient", "Disease")], fractions)

saveRDS(df, "results/processed/plasma_cell_fractions_idol_with_metadata.rds")

library(dplyr)

df_summary <- df %>%
  group_by(Disease) %>%
  summarise(across(Bcell:NK, list(mean = mean, sd = sd)))

print(df_summary)

library(tidyr)

df_mean <- df %>%
  group_by(Disease) %>%
  summarise(across(Bcell:NK, mean)) %>%
  pivot_longer(
    cols = Bcell:NK,
    names_to = "CellType",
    values_to = "Fraction"
  )

library(ggplot2)

print(ggplot(df_mean, aes(x = Disease, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "Average immune composition by disease",
    y = "Mean fraction",
    fill = "Cell type"
  )
)

print(kruskal.test(Fraction ~ Disease, data = subset(df_mean, CellType == "Mono")))

library(pheatmap)

mat <- as.matrix(df[, 4:10])
rownames(mat) <- df$sample_name

ann_colors <- list(
  Disease = c(
    MG  = "#1f77b4",
    SMM = "#ff7f0e",
    MM  = "#d62728"
  )
)


pheatmap(
  fractions,
  scale = "row",
  annotation_row = df["Disease", drop = FALSE],
  annotation_colors = ann_colors,
  clustering_method = "ward.D2",
  main = "Plasma immune composition"
)

library(ggfortify)

pca <- prcomp(df[, 4:10], scale. = TRUE)

autoplot(
  pca,
  data = df,
  colour = "Disease",
  size = 3
)

cor_mat <- cor(df[, 4:10])
pheatmap(cor_mat, main = "Cell-type correlation matrix")

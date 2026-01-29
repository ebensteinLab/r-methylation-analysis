#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

library(tidyverse)

BASE_DIR <- "raw_data/cellfree"
SAMPLESHEET <- file.path(BASE_DIR, "samplesheet.csv")

# -----------------------------
# Load global sample sheet
# -----------------------------
if (!file.exists(SAMPLESHEET)) {
  stop("Global samplesheet.csv not found in ", BASE_DIR)
}

cf_targets <- read.csv(
  SAMPLESHEET,
  stringsAsFactors = FALSE,
  colClasses = "character"
)

message("Loaded ", nrow(cf_targets), " samples from global sample sheet")

# -----------------------------
# Construct Directory and Basename fields
# -----------------------------
cf_targets$Directory <- file.path(BASE_DIR, cf_targets$Sentrix_Id)

cf_targets$Basename <- file.path(
  cf_targets$Directory,
  paste0(cf_targets$Sentrix_Id, "_", cf_targets$Sentrix_Position)
)

# -----------------------------
# Validate IDAT existence
# -----------------------------
cf_targets$GrnExists <- file.exists(paste0(cf_targets$Basename, "_Grn.idat"))
cf_targets$RedExists <- file.exists(paste0(cf_targets$Basename, "_Red.idat"))

missing <- cf_targets %>% filter(!GrnExists | !RedExists)

if (nrow(missing) > 0) {
  message("❌ Missing IDAT files detected:")
  print(missing[, c("Sentrix_Id", "Sentrix_Position", "Basename")])
  stop("Fix missing IDAT files before proceeding.")
}

message("All IDAT files found ✔")

# Drop helper columns
cf_targets <- cf_targets %>% select(-GrnExists, -RedExists)

# -----------------------------
# Save merged cf_targets
# -----------------------------
dir.create("results/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(cf_targets, "results/processed/cf_targets_merged.rds")

rm(cf_targets)
gc()

message("Saved cf_targets to results/processed/cf_targets_merged.rds")

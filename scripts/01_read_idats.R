#!/usr/bin/env Rscript

library(tidyverse)

# -----------------------------
# Directories containing sample sheets
# -----------------------------
dirs <- c(
  "raw_data/209547710016", 
  "raw_data/209742170040", 
  "raw_data/209742170042", 
  "raw_data/209742180031", 
  "raw_data/209742180157"
)

# -----------------------------
# Load and merge sample sheets
# -----------------------------
targets <- bind_rows(lapply(dirs, function(d) {
  ss <- read.csv(
    file.path(d, "samplesheet.csv"),
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  ss$Directory <- d
  ss
}))

# -----------------------------
# Construct Basename for SeSAMe
# -----------------------------
targets$Basename <- file.path(
  targets$Directory,
  paste0(targets$Sentrix_Id, "_", targets$Sentrix_Position)
)

# -----------------------------
# Validate existence of IDAT files
# -----------------------------
targets$GrnExists <- file.exists(paste0(targets$Basename, "_Grn.idat"))
targets$RedExists <- file.exists(paste0(targets$Basename, "_Red.idat"))

if (any(!targets$GrnExists | !targets$RedExists)) {
  message("WARNING: Some IDAT files are missing!")
  print(targets %>% filter(!GrnExists | !RedExists))
  stop("Fix missing IDAT files before proceeding.")
}

# Remove those helper columns before saving
targets <- targets %>% select(-GrnExists, -RedExists)

# -----------------------------
# Save merged sample sheet
# -----------------------------
saveRDS(targets, "results/processed/targets_merged.rds")
message("Saved: results/processed/targets_merged.rds")


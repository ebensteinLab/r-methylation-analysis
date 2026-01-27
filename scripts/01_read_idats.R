#!/usr/bin/env Rscript

if (!endsWith(getwd(), "R/projects/r-methylation-analysis")) {
  setwd("R/projects/r-methylation-analysis")
}

library(tidyverse)

# -----------------------------
# Directories containing sample sheets
# -----------------------------
find_idat_dirs <- function(
    base_dir = "raw_data/genomic",
    samplesheet = "samplesheet.csv",
    recursive = FALSE
) {
  # List candidate directories
  dirs <- list.dirs(base_dir, recursive = recursive, full.names = TRUE)
  dirs <- dirs[dirs != base_dir]
  
  valid_dirs <- character()
  
  for (d in dirs) {
    idats <- list.files(d, pattern = "\\.idat?$", ignore.case = TRUE)
    
    # Skip dirs without IDATs
    if (length(idats) == 0) {
      next
    }
    
    sheet_path <- file.path(d, samplesheet)
    
    if (!file.exists(sheet_path)) {
      warning(
        sprintf(
          "Directory '%s' contains IDAT files but is missing %s — skipping.",
          d, samplesheet
        ),
        call. = FALSE
      )
      next
    }
    
    valid_dirs <- c(valid_dirs, d)
  }
  
  message("Found ", length(valid_dirs), " valid IDAT directories.")
  return(valid_dirs)
}

dirs <- find_idat_dirs()

print(dirs)

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

rm(targets)
gc()
message("Saved: results/processed/targets_merged.rds")


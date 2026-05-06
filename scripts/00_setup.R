#!/usr/bin/env Rscript

# ----------------------------------------------------------------------
# Script: 00_setup.R (Updated for Stable R 4.4 / Bioc 3.20)
# ----------------------------------------------------------------------

message("Starting setup for Stable R 4.4 / Bioc 3.20...")

# -----------------------------
# 1. Set System Environment & Timeouts
# -----------------------------
# Increase timeout for large annotation and data downloads
options(timeout = 1000) 
options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

# -----------------------------
# 2. Pin Bioconductor Version
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Force installation of Bioconductor 3.20 (Stable)
# This is critical to avoid the metadata errors seen in 3.21/3.22
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
  
# -----------------------------
# 3. Helper Functions
# -----------------------------
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_bioc_if_missing <- function(pkgs) {
  options(repos = BiocManager::repositories())
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE, version = "3.20")
    }
  }
}

# -----------------------------
# 4. Install Core CRAN Packages
# -----------------------------
cran_pkgs <- c("tidyverse", "data.table", "matrixStats", "ggplot2", 
               "RColorBrewer", "nnls", "compositions", "remotes")
install_if_missing(cran_pkgs)

# -----------------------------
# 5. Install Core Bioconductor Packages
# -----------------------------
bioc_pkgs <- c("ExperimentHub", "AnnotationHub", "BiocParallel", 
               "sva", "limma", "EpiDISH", "sesame", "sesameData",
               "IlluminaHumanMethylationEPICmanifest",
               "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
               "GenomicRanges", "DMRcate", "GEOquery")
install_bioc_if_missing(bioc_pkgs)

# -----------------------------
# 7. Initialize Caches
# -----------------------------
# Pre-create and clear any old/corrupt development caches
eh_cache <- path.expand("~/.cache/R/ExperimentHub")
dir.create(eh_cache, recursive = TRUE, showWarnings = FALSE)

# Cache SeSAMe data
if (requireNamespace("sesameData", quietly = TRUE)) {
  sesameData::sesameDataCacheAll()
}

# -----------------------------
# 8. Directory Structure
# -----------------------------
dirs <- c("results/processed", "results/deconvolution", "results/qc", "scripts")
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 6. Install Deconvolution Reference
# -----------------------------
# Install from GitHub to ensure we have the latest EPICv2-compatible code
# but within the Bioc 3.20 environment
remotes::install_github("immunomethylomics/FlowSorted.BloodExtended.EPIC", 
                        dependencies = TRUE, 
                        upgrade = "never")

message("\nSetup complete. Environment is now pinned to Bioconductor 3.20.")


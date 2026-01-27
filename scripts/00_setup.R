#!/usr/bin/env Rscript

# ----------------------------------------------------------------------
# Script: 00_setup.R
# Purpose: Cleanly install all required libraries for SeSAMe EPICv2 processing
# ----------------------------------------------------------------------

# -----------------------------
# System prerequisites (LinuxgetOption("repos"))
# -----------------------------
# sudo apt update
# sudo apt install -y libharfbuzz-dev libfribidi-dev libfreetype6-dev \
#     libpng-dev libfontconfig1-dev libcurl4-openssl-dev libxml2-dev libssl-dev

message("Starting setup...")

# -----------------------------
# Helper: Install CRAN packages
# -----------------------------
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing CRAN package: %s", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
      message(sprintf("Package already installed: %s", pkg))
    }
  }
}

# -----------------------------
# Helper: Install Bioconductor packages
# -----------------------------
install_bioc_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  options(repos = BiocManager::repositories())
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing Bioconductor package: %s", pkg))
      BiocManager::install(pkg, ask = FALSE, update = TRUE)
    } else {
      message(sprintf("Bioc package already installed: %s", pkg))
    }
  }
}

# -----------------------------
# **IMPORTANT** Remove old/broken SeSAMe installations
# (These cause the missing-guessPlatform issue.)
# -----------------------------
#remove_old_sesame <- function() {
#  lib_paths <- .libPaths()
#  for (lib in lib_paths) {
#    sesame_path <- file.path(lib, "sesame")
#    sesameData_path <- file.path(lib, "sesameData")
#    if (dir.exists(sesame_path)) {
#      message(sprintf("Removing old sesame installation at: %s", sesame_path))
#      unlink(sesame_path, recursive = TRUE, force = TRUE)
#    }
#    if (dir.exists(sesameData_path)) {
#      message(sprintf("Removing old sesameData installation at: %s", sesameData_path))
#      unlink(sesameData_path, recursive = TRUE, force = TRUE)
#    }
#  }
#}

#remove_old_sesame()

# -----------------------------
# Required CRAN packages
# -----------------------------
cran_pkgs <- c(
  "tidyverse",
  "data.table",
  "matrixStats",
  "ggplot2",
  "RColorBrewer",
  "remotes"     # used for GitHub fallback installs
)

install_if_missing(cran_pkgs)

remotes::install_github("qsbase/qs")

# -----------------------------
# Install SeSAMe and sesameData cleanly
# -----------------------------
bioc_pkgs <- c(
  "ExperimentHub",
  "AnnotationHub",
  "BiocParallel",
  "sva",
  "limma"
)

install_bioc_if_missing(bioc_pkgs)

plot_pkgs <- c(
  "pheatmap",
  "ggrotify"
)

install_bioc_if_missing(plot_pkgs)

# -----------------------------
# Step 1: Try installing official Bioconductor SeSAMe
# -----------------------------
message("Installing SeSAMe from Bioconductor...")
try({
  BiocManager::install("sesame", ask = FALSE, update = TRUE)
  BiocManager::install("sesameData", ask = FALSE, update = TRUE)
}, silent = TRUE)


library(sesame)

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
  "GenomicRanges",
  "DMRcate"
))

# -----------------------------
# Avoid interactive ExperimentHub/AnnotationHub prompts
# -----------------------------
options(ExperimentHub.ask = FALSE)
options(AnnotationHub.ask = FALSE)

# Pre-create cache directories
eh <- Sys.getenv("EXPERIMENT_HUB_CACHE", unset = "~/.cache/R/ExperimentHub")
eh <- path.expand(eh)
ah <- "~/.cache/R/AnnotationHub"
dir.create(eh, recursive = TRUE, showWarnings = FALSE)
dir.create(ah, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Recommended extras
# -----------------------------
optional_bioc <- c(
  "GEOquery",
  "ewastools"
)

install_bioc_if_missing(optional_bioc)

sesameDataCacheAll()

# -----------------------------
# Project directory structure
# -----------------------------
dirs <- c(
  "results",
  "results/processed",
  "results/features",
  "results/qc",
  "scripts"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("Created directory: %s", d))
  }
}

# -----------------------------
# Session Info
# -----------------------------
message("\nAll setup complete. Session Info:")
print(sessionInfo())

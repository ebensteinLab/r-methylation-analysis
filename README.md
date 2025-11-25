# r-methylation-analysis
R analysis pipeline of methylation data

install.packages("BiocManager")
BiocManager::install(c(
  "minfi",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "GEOquery",
  "sva",
  "limma"
))

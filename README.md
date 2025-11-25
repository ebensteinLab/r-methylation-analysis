# r-methylation-analysis
R analysis pipeline of methylation data

install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install(c("rtracklayer", "biomaRt"))
BiocManager::install("GenomicFeatures")
BiocManager::install("sva")
BiocManager::install(c("minfi","IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))

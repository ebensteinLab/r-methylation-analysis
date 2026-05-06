library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Load annotation
ann <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# ann has columns: chr, pos
epic_df <- data.frame(
  probe_id = rownames(ann),
  chrom = ann$chr,
  pos = ann$pos,
  stringsAsFactors = FALSE
)

# Load your lifted CpGs (hg38)
lifted <- read.delim("sandbox/singles_blood.tsv")

# Merge
matches <- merge(
  lifted,
  epic_df,
  by.x = c("chrom", "central_pos"),
  by.y = c("chrom", "pos"),
  all.x = TRUE
)

write.csv(matches, "sandbox/singles_blood_hg38.csv")

rm(ann, epic_df, matches, lifted)

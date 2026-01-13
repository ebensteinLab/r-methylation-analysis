# r-methylation-analysis
R analysis pipeline of methylation data

Before openning RStudio, create file ~/.Rprofile, and paste the following lines in it:
Sys.setenv(
  EXPERIMENT_HUB_CACHE = path.expand("~/.cache/sesameData"),
  SESAMEDATA_CACHE = path.expand("~/.cache/sesameData")
)
If you already opened RStudio before adding the above lines to the file, you
need to close RStudio, add the lines to the file and open RStudio again.

Go to the project root directory in command line, and perform these two commands:
- wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
- gunzip hg19ToHg38.over.chain.gz

Run in R studio the scripts in the "scripts" directory according to their numbers

Requires R version >= 4.5


library(ShortRead)

base <- "FKR02_3bfmi_R2_001"

f1 <- FastqStreamer(paste0(base, ".fastq.gz"), n=25000000)
batch <- 1
while (length(fq1 <- yield(f1))) {
  subsetname <- paste0("subset/", base, ".subset", batch, ".fastq.gz")
  cat("Split into:", subsetname, "\n")
  writeFastq(fq1, subsetname)
  batch <- batch + 1
}
close(f1)

## RNA Seq analysis pipeline, subset fastq files to only contain 500.000 reads
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2015
# first written Jun, 2015

### Read the sample description file

library(ShortRead)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

samples <- read.table("Experiment/SampleDescription.txt", sep="\t", header=TRUE)
sample.names <- apply(samples[, c("Lib_id","TagLane")], 1, paste0, collapse="_")
names(sample.names) <- samples[, "Lib_id"]

for(name in sample.names){
  if(file.exists(paste0("FASTQ/", name, "_R1_001.fastq.gz"))){
    cat("Creating subset", name, "\n")
    f1 <- FastqSampler(paste0("FASTQ/", name, "_R1_001.fastq.gz"), n=500000)
    f2 <- FastqSampler(paste0("FASTQ/", name, "_R2_001.fastq.gz"), n=500000)
    set.seed(123L); p1 <- yield(f1);
    set.seed(123L); p2 <- yield(f2);
    writeFastq(p1, paste0("FASTQ/", name, "_R1_001.subset.fastq.gz"))
    writeFastq(p2, paste0("FASTQ/", name, "_R2_001.subset.fastq.gz"))
  }else{
    cat("Skipping", name, "\n")
  }
}

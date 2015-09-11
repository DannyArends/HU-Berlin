# Generate commands for the RNAseq pipeline, for all the samples in sampledescription.txt
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sept, 2015
# first written Juli, 2015

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/Experiment")
samples <- read.table("sampledescription.txt", sep="\t", header=TRUE, colClasses="character")

for(x in 1:nrow(samples)){
  fname <- paste0(samples[x,"Lib_id"], "_", samples[x,"TagLane"])
  cat(paste0("nohup Rscript pipeline.R /home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/FASTQ/", fname," > ",fname,".log 2>&1&"),"\n")
}

# Re-analysis of the micro array data from high fat and low fat chickens
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

setwd("E:/Chicken/RNA/Arrays_for_fat-2011-01/")
arrays <- read.table("Annotation/arrays.txt", sep="\t", header=TRUE, colClasses = "character")
lowfat  <- arrays[which(arrays[,"fat"] == "low"),"filename"]
highfat <- arrays[which(arrays[,"fat"] == "high"),"filename"]

group1 <- paste0("lowfat",1:4)
group2 <- paste0("highfat",1:4)

fatdata <- NULL
for(filename in c(lowfat,highfat)){
  mdata <- read.csv(filename, sep = "\t", skip = 9, header = TRUE)
  fatdata <- cbind(fatdata, mdata$gProcessedSignal)
}

fatdata <- cbind(mdata[,c("ProbeName", "Sequence")], fatdata)
colnames(fatdata) <- c("ProbeName", "Sequence", group1, group2)

### Create fasta file

# Code see CJ

# blastn -task blastn -query E:\Chicken\RNA\Arrays_for_fat-2011-01\Analysis\probes.fasta -db Gallus_gallus.Galgal4.74.dna.db -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 5 -out E:\Chicken\RNA\Arrays_for_fat-2011-01\Analysis\ProbeLocations.txt


### Quality control !!!

### Create a fasta file, and blast to the genome, re annotate the data that we have

### Quantile normalisation of the data

### Find differentially expressed probes

### Summarize this into genes
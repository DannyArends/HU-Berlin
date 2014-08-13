# SNPanalysisDiabetes.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written May, 2014
#
# Analysis of data from Sebastiaan

setwd("E:/Mouse/DNA/DiversityArray/")
SNPdata <- read.table("Analysis/measurementsAtlas_annotated.txt", sep="\t", header=TRUE, colClasses="character")

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                           # Load the annotation
validAnnot <- match(colnames(SNPdata)[9:ncol(SNPdata)], rownames(annotation))                     # Annotation that matches the mice we have in our SNP data
annotation <- annotation[validAnnot, ]

header <- cbind(matrix("",7,8),t(annotation))                                                     # Create the table header
colnames(header)[1:8] <- colnames(SNPdata)[1:8]                                                   # Add column names to the empty columns

chromosomes <- c(1:19, "X", "Y", "MT")
for(chr in chromosomes){
  write.table(rbind(header, SNPdata[which(SNPdata[,"Chr"] == chr),]), paste0("Analysis/perChromosome/MouseDiversity_",chr,".txt"), sep="\t", row.names = FALSE)
}

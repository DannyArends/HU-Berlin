# Analysis of the micro array data from Atlas 2014
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

region <- c(3, 36474653, 36854743)  # Determined from the MegaMuga: JAX00106222 to JAX00519845

setwd("E:/Mouse/RNA/ArrayDesign/Atlas data")
arrays <- read.table("Annotation/arrays.txt", header=TRUE, sep="\t", colClasses="character")
alldata <- read.table("Analysis/geneexpression.txt", sep="\t", header=TRUE)

onchr        <- which(alldata[,"chromosome_name"] == region[1])
insideregion <- which(alldata[,"chromosome_name"] == region[1] & alldata[,"start_position"] > region[2] & alldata[,"end_position"] < region[3])

cat("Found", length(onchr),"probes that target", length(unique(alldata[onchr,"ensembl_gene_id"])),"unique genes on our chromosome of interest\n")
cat("Found", length(insideregion),"probes that target", length(unique(alldata[insideregion,"ensembl_gene_id"])),"unique genes inside our region of interest\n")

write.table(alldata[onchr,-11], "Analysis/genesOnChromosome3.txt", sep="\t", row.names=FALSE)
write.table(alldata[insideregion,-11], "Analysis/genesInRegion.txt", sep="\t", row.names=FALSE)

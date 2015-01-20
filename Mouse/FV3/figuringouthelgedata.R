# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

library(biomaRt)                                                                                                      # Biomart - Used to annotate probes
library(preprocessCore)                                                                                               # preprocessCore - For quantile normalisation

setwd("E:/Mouse/RNA/FV3")

arrays    <- read.table("Annotation/arrays.txt", sep="\t", header=TRUE, row.names=1)                                   # Arrays annotation
probes    <- read.table("Annotation/GPL6105-11577.txt", sep="\t", header=TRUE, row.names=1)                                   # Arrays annotation

rawdata   <- read.table("RawData/expressions.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)              # Load the expression data
normdata  <- read.table("Analysis/NormData.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)              # Load the expression data


which(!rownames(rawdata) %in% probes[,"Array_Address_Id"])
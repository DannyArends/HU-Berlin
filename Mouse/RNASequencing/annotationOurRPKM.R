# annotations.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014
#
# Annotating the RNA sequencing data done on our own server

library(biomaRt)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

RPKM <- read.csv("Analysis/RPKM.txt", sep="\t", header=TRUE, check.names=FALSE)                                                          # RNA-Seq RPKM data
annotation <- read.csv("FastQ/sampledescription.txt", sep="\t", header=TRUE)                                             # Sample annotation
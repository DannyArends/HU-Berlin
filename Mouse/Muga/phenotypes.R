# Muga Phenotypes
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/DNA/MegaMuga/")                                                             # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t")
genotypes <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                              # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)


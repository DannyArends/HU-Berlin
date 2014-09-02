# Muga Genotypes
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/DNA/MegaMuga/")                                                             # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)

numgeno <- apply(genotypes, 2, function(x){as.numeric(as.factor(x))})

dendrogram <- as.dendrogram(hclust(dist(t(numgeno))))

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
missingPerMar <- apply(genotypes, 1, function(x){ sum(is.na(x)) / length(x) * 100 })

genotypes <- genotypes[,-which(missingPerInd==100)]
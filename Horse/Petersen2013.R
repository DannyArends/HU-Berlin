# Analysis of horse SNP chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2016
# first written Jan, 2016

setwd("E:/Horse/DNA/Petersen2013")

description <- read.table("KSF1385391590.txt", header = FALSE, sep = "\t")
colnames(description) <- c("shortname", "fullname")

map <- read.table("AZK1385391535.map", row.names = 2)[,-2]

genotypedata <- read.table("JXB1385391482.ped", na.strings = c("NA", "-9", "0"))
individuals <- genotypedata[,1:6]
genotypedata <- genotypedata[,7:ncol(genotypedata)]

coding <- c("A","C","T","G")
x <- 1; paste0(coding[genotypedata[,x]],coding[genotypedata[,x+1]])

genotypes <- matrix(NA, dim(genotypedata)[1], nrow(map))
cnt <- 1
for(x in seq(1, ncol(genotypedata),2)){
  na <- c(which(is.na(genotypedata[,x])),which(is.na(genotypedata[,x+1])))
  marker <- paste0(coding[genotypedata[,x]],coding[genotypedata[,x+1]])
  if(length(na >= 1)) marker[na] <- NA
  genotypes[,cnt] <- marker
  cnt <- cnt + 1
}

rownames(genotypes) <- paste0(individuals[,1],"_", individuals[,2])
colnames(genotypes) <- rownames(map)

# Load in our data

setwd("E:/Horse/DNA/Equine60k/")
mapA <- read.table(file="input/cleaned_map.txt", sep = "\t")
genotypesA <- read.table(file="input/cleaned_genotypes.txt", sep = "\t")                          # Save the clean genotypes to disk

shared <- which(rownames(map) %in% paste0("chr", mapA[,1],".",mapA[,2]))

map <- map[shared,]
genotypes <- genotypes[,shared]

shared <- which(paste0("chr", mapA[,1],".",mapA[,2]) %in% rownames(map))

mapA <- mapA[shared,]
genotypesA <- genotypesA[shared,]

rownames(mapA) <- paste0("chr", mapA[,1],".",mapA[,2])
rownames(genotypesA) <- paste0("chr", mapA[,1],".",mapA[,2])

# Combine the data
mapA <- mapA[match(rownames(map), rownames(mapA)),]
genotypesA <- genotypesA[match(rownames(map), rownames(genotypesA)),]






arabianHorses <- genotypes[which(grepl("_ARR", rownames(genotypes))),]

# Merge ASE data from 3 tissues
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross/Analysis")
Gf <- read.csv("snp_stats_Gonadal Fat.txt", sep = "\t", colClasses="character")
rownames(Gf) <- apply(Gf[,c(1:4,7)], 1,paste0, collapse="_")
colnames(Gf)[13:ncol(Gf)] <- paste0(colnames(Gf)[13:ncol(Gf)], "_GonadalFat")

Li <- read.csv("snp_stats_Liver.txt", sep = "\t", colClasses="character")
rownames(Li) <- apply(Li[,c(1:4,7)], 1,paste0, collapse="_")
colnames(Li)[13:ncol(Li)] <- paste0(colnames(Li)[13:ncol(Li)], "_Liver")

Qa <- read.csv("snp_stats_Quadriceps.txt", sep = "\t", colClasses="character")
rownames(Qa) <- apply(Qa[,c(1:4,7)], 1,paste0, collapse="_")
colnames(Qa)[13:ncol(Qa)] <- paste0(colnames(Qa)[13:ncol(Qa)], "_Quadriceps")

uniqueRows <- unique(c(rownames(Gf), rownames(Li), rownames(Qa)))
uniqueColumns <- unique(c(colnames(Gf), colnames(Li), colnames(Qa)))

combined <- data.frame(matrix(NA, length(uniqueRows), length(uniqueColumns), dimnames=list(uniqueRows, uniqueColumns)))

for(x in rownames(combined)){
  if(x %in% rownames(Gf)) combined[x, colnames(Gf)] <- unlist(Gf[x, ])
  if(x %in% rownames(Li)) combined[x, colnames(Li)] <- unlist(Li[x, ])
  if(x %in% rownames(Qa)) combined[x, colnames(Qa)] <- unlist(Qa[x, ])
}

write.table(combined, "combined_ASE_3tissues.txt", sep = '\t')
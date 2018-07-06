#
# Analysis of BxD genotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

setwd("E:/Mouse/BxD/CTL")

genotypes <- read.table("BXD.geno",sep="\t", skip=6, header=TRUE, row.names=2, na.strings="U")
map <- genotypes[,1:3]
genotypes <- genotypes[,-c(1:3)]

map[1:10,]
genotypes[1:10,]

GN379 <- read.csv("GN379_MeanDataAnnotated_rev081815.txt", skip=33, header=TRUE, sep="\t")
GN379[1:10,]

renames <- c("BXD96","BXD97","BXD92","BXD80","BXD103")
names(renames) <- c("BXD48a", "BXD65a", "BXD65b", "BXD73a", "BXD73b")

for(i in which(colnames(GN379) %in% names(renames))){
  colnames(GN379)[i] <- renames[colnames(GN379)[i]]
}

measured <- colnames(GN379)[which(grepl("BXD", colnames(GN379)))]
genotypes <- genotypes[,measured]

numgeno <- apply(genotypes, 1, function(x){
  as.numeric(as.factor(x))
})

numgeno[is.na(numgeno)] <- -9
GN379 <- t(GN379[,measured])

library(ctl)

s <- proc.time()
res <- CTLscan(numgeno, GN379, phenocol = 1, verbose = TRUE)
(proc.time() - s)[3]


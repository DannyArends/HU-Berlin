#
# Smoothing genotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

library("openxlsx")
source("D:/Github/HU-Berlin/RobWilliams/functions.R")

setwd("E:/Mouse/BxD")
alldata <- read.csv("genotypes.txt", sep = "\t", colClasses="character", skip = 1, row.names=2, header=TRUE, na.strings=c("NA", ""))

map <- alldata[,c(3, 5, 7)]                                                       # Extract the physical map
genotypes <- alldata[,-c(1:11)][,1:203]                                           # Extract the genotypes

genotypes[3238, "BXD53"]  <- NA
genotypes[3239, "BXD53"]  <- NA
genotypes[3240, "BXD53"]  <- NA
genotypes[3241, "BXD53"]  <- NA
genotypes[20364,"BXD146"] <- NA
genotypes[16041,"BXD93"]  <- NA
genotypes[16042,"BXD93"]  <- NA

#table(unlist(genotypes))

toNumGeno <- function(genotypes){
  numgeno <- apply(genotypes, 2, function(x){ return(as.numeric(as.character(x))) })      # transform genotypes to numeric values
  rownames(numgeno) <- rownames(genotypes)                                                # Set the markernames
  return(numgeno)
}

geno <- toNumGeno(genotypes)

geno <- geno[-which(rownames(geno) == "Affy_10542764_Arntl2"),]             # Marker without a locations
map <- map[-which(rownames(map) == "Affy_10542764_Arntl2"),]                # Marker without a locations

#table(unlist(geno))

recomb.index <- apply(geno, 2, function(x){recombinations(x, chr = map[,1]) })
recomb.loc <- apply(geno, 2, function(x){recombinations(x, loc = as.numeric(map[,2]), chr = map[,1], TRUE) })

wb <- createWorkbook()
sheet <- addWorksheet(wb, sheetName="Genotypes BXD")
writeData(wb, sheet = 1, rownames(geno), startCol = 1, startRow = 2)
writeData(wb, sheet = 1, map[rownames(geno),1], startCol = 2, startRow = 2)
writeData(wb, sheet = 1, map[rownames(geno),2], startCol = 3, startRow = 2)
writeData(wb, sheet = 1, map[rownames(geno),3], startCol = 4, startRow = 2)
for(x in 1:ncol(geno)){
  writeData(wb, sheet = 1, c(colnames(geno)[x], geno[,x]), startCol = (x+4))
}

recombOK <- createStyle(fontColour = rgb(1,1,1), fgFill = rgb(0,0,1))
recombBAD <- createStyle(fontColour = rgb(1,1,1), fgFill = rgb(1,0,0))

for(i in 1:ncol(geno)){
  rows <- as.numeric(unlist(apply(cbind(floor(as.numeric(recomb.index[[i]])), ceiling(as.numeric(recomb.index[[i]]))), 1, function(x){ return(x[1]:x[2])})))
  addStyle(wb, sheet = 1, recombOK, cols = rep((i+4), length(rows)), rows = (rows+1))
  diffs <- diff(as.numeric(get.locs(recomb.loc[[i]])))
  chrs  <- diff(as.numeric(as.factor(get.chr(recomb.loc[[i]]))))
  BADlocsFrom <- recomb.index[[i]][which(diffs > 0 & diffs < 1 & chrs == 0)]
  BADlocsTo <- recomb.index[[i]][which(diffs > 0 & diffs < 1 & chrs == 0) + 1]
  rows <- as.numeric(unlist(apply(cbind(floor(as.numeric(BADlocsFrom)), ceiling(as.numeric(BADlocsTo))), 1, function(x){ return(x[1]:x[2])})))
  addStyle(wb, sheet = 1, recombBAD, cols = rep((i+4), length(rows)), rows = (rows+1))
}

saveWorkbook(wb, file="genotypes_recombinations.xlsx", overwrite = TRUE)

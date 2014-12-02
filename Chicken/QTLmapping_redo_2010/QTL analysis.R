# QTL analysis of the QTL data from Mustapha
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

setwd("E:/Chicken/ClassicalPhenotypes/QTLanalysis_Mostafa_2010")

genotypes <- read.table("genotypes.txt", sep="\t", na.strings = ".", header=TRUE)
bodycomposition <- read.table("bodycomposition.txt", sep="\t", na.strings = ".", header=TRUE)[, -1]

BCinGT <- which(bodycomposition[,"ID"] %in% genotypes[,"ID"])
bodycomposition <- bodycomposition[BCinGT, ]

GTinBC <- which(genotypes[,"ID"] %in% bodycomposition[,"ID"])
genotypes <- genotypes[GTinBC, ]

GTorderToMatchBC <- match(as.character(bodycomposition[,"ID"]), as.character(genotypes[,"ID"]))
genotypes <- genotypes[GTorderToMatchBC, ]

rownames(genotypes) <- genotypes[,"ID"]
rownames(bodycomposition) <- bodycomposition[,"ID"]

mapQTL <- function(bodycomposition, genotypes, phenotypename = "Body.weight.at.20.weeks..g"){
  markersForQTL <- colnames(genotypes)[7:ncol(genotypes)]
  LODscores <- NULL
  for(markername in markersForQTL){
    res <- anova(lm(bodycomposition[,phenotypename] ~  bodycomposition[,"HatchDate"] + bodycomposition[,"Sire"] + genotypes[,markername]))
    LODscores <- rbind(LODscores, -log10(res[[5]]))
  }
  rownames(LODscores) <- markersForQTL
  colnames(LODscores) <- c("HatchDate", "Subfamily", "Marker", "Residuals")
  return(LODscores)
}

phenotypesForQTL <- colnames(bodycomposition)[8:ncol(bodycomposition)]

LODresults <- vector("list", length(phenotypesForQTL))
names(LODresults) <- phenotypesForQTL
for(phenotype in phenotypesForQTL){
  LODresults[[phenotype]] <- mapQTL(bodycomposition, genotypes, phenotype)
}

for(x in 1:length(LODresults)){
  above4 <- which(LODresults[[x]][,"Marker"] > 4.635)
  cat(names(LODresults)[x], names(above4),"\n")
}

## What do we learn: the first 6 traits show a similar QTL profile from:
## Small region: rs14488203 => 4 69583407 to rs14494169 => 4 78715886
## Extended region: MCW0276 => 4 59553432 to UMA4024 =>    4 84774762

plot(LODresults$Mesenteric.adipose.tissue..g[,"Marker"], t='l')

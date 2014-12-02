# QTL analysis of the QTL data from Mustapha
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

setwd("E:/Chicken/ClassicalPhenotypes/QTLanalysis_Mostafa_2010")

genotypes <- read.table("genotypes.txt", sep="\t", na.strings = ".", header=TRUE)                                 # Load the genotypes
bodycomposition <- read.table("bodycomposition.txt", sep="\t", na.strings = ".", header=TRUE)[, -1]               # Load the phenotypes

BCinGT <- which(bodycomposition[,"ID"] %in% genotypes[,"ID"])                                                     # Which body composition traits have genotypes
bodycomposition <- bodycomposition[BCinGT, ]

GTinBC <- which(genotypes[,"ID"] %in% bodycomposition[,"ID"])                                                     # Which genotypes have body composition traits
genotypes <- genotypes[GTinBC, ]

GTorderToMatchBC <- match(as.character(bodycomposition[,"ID"]), as.character(genotypes[,"ID"]))                   # Reorder the genotypes such that they match the body composition ordering
genotypes <- genotypes[GTorderToMatchBC, ]

rownames(genotypes) <- genotypes[,"ID"]
rownames(bodycomposition) <- bodycomposition[,"ID"]

mapQTL <- function(bodycomposition, genotypes, phenotypename = "Body.weight.at.20.weeks..g"){                     # Function to do QTL mapping using 2 covariates
  markersForQTL <- colnames(genotypes)[7:ncol(genotypes)]                                                         # Markers that are used in QTL mapping
  LODscores <- NULL                                                                                               # Matrix for results
  for(markername in markersForQTL){
    res <- anova(lm(bodycomposition[,phenotypename] ~  bodycomposition[,"HatchDate"] + bodycomposition[,"Sire"] + genotypes[,markername]))
    LODscores <- rbind(LODscores, -log10(res[[5]]))
  }
  rownames(LODscores) <- markersForQTL
  colnames(LODscores) <- c("HatchDate", "Subfamily", "Marker", "Residuals")                                      # Add the row names to our result matrix
  return(LODscores)
}

phenotypesForQTL <- colnames(bodycomposition)[8:ncol(bodycomposition)]                                           # Phenotypes that are used in QTL mapping

LODresults <- vector("list", length(phenotypesForQTL))                                                           # Empty list of QTL matrices, one for each phenotype
names(LODresults) <- phenotypesForQTL
for(phenotype in phenotypesForQTL){
  LODresults[[phenotype]] <- mapQTL(bodycomposition, genotypes, phenotype)                                       # Map the QTL and store the resulting matrix in the list
}

for(x in 1:length(LODresults)){
  above4 <- which(LODresults[[x]][,"Marker"] > 4.65)                                                             # Analyse the resulting profiles, and look at which markers are above the threshold
  cat(names(LODresults)[x], names(above4),"\n")                                                                  # -log10(0.05/ (123 * 18)) = 4.65
}
plot(LODresults$Mesenteric.adipose.tissue..g[,"Marker"], t='l')                                                  # Plot one of our QTL profiles

## What do we learn: the first 6 traits show a similar QTL profile from:
## Small region: rs14488203 => 4 69583407 to rs14494169 => 4 78715886
## Extended region: MCW0276 => 4 59553432 to UMA4024 =>    4 84774762

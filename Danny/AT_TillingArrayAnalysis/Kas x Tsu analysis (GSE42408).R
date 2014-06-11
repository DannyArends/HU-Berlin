# Kas x Tsu analysis (GSE42408).R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2014
# first written May, 2014
#
# Script for pre-processing, log transformation, normalisation, and CTL/QTL analysis of Kas x Tsu data (GSE42408)

# SERVER PART === THIS NEEDS A LOT OF RAM
setwd("/home/arends")
library(affy)
affy <- ReadAffy(filenames=paste0("Data/", dir("Data")))

nArrayCols <- 1612

# Sequences & mappings of the probes, they seem to be mapped against TAIR 10
# probeSequences <- read.table("Meta/GPL16303_TilingatSNPtilx520433_At_TAIRG_probe_tab.txt", sep="\t", header=TRUE)

probeMappings <- read.table("Meta/GPL16303_TilingatSNPtilx520433_At_TAIRG_mapping.txt", sep="\t", header=TRUE)

indices <- NULL                                           # Translate the x,y coordinates to array indices
for(r in 1:nrow(probeMappings)){
  indices <- c(indices, xy2indices(x = probeMappings[r,"Probe.X"], y = probeMappings[r,"Probe.Y"], nc = nArrayCols))
  if(r %% 1000 == 0){
    cat(paste0("Done ",r,"\n")); gc();
  }
}

allexpressions <- exprs(affy)
probeexpressions <- allexpressions[indices, ]
annotatedProbes <- cbind(probeMappings[, 1:4], probeexpressions)
write.table(annotatedProbes, file="annotatedExpressionData.txt", sep = "\t", row.names = FALSE)

# END OF SERVER PART
# LOCAL

setwd("d:/CTLdata")
library(limma)

rawdata <- read.table("Kas_Tsu_RIdata.txt", sep='\t', header=TRUE, check.names=FALSE)             # Load our genotypes
envdata <- read.table("environments.txt", sep=' ', header=FALSE, row.names=1)                     # Environments

mapraw              <- t(rawdata[1:2, 2:ncol(rawdata)])                                           # Store the genetic map
map                 <- apply(mapraw,2,as.numeric)                                                 # Numeric distances
rownames(map)       <- rownames(mapraw)                                                           # Names of the markers

genotypes           <- rawdata[3:nrow(rawdata), 2:ncol(rawdata)]                                  # Remove the map
rownames(genotypes) <- rawdata[3:nrow(rawdata),1]                                                 # Names of the individuals

genotypes <- genotypes[-c(1, 2),]                                                                 # Remove the founders

genonum <- apply(genotypes, 2, function(x){ as.numeric(as.factor(x)) })                           # Numeric genotypes
rownames(genonum) <- rownames(genotypes)                                                          # Names of the individuals

individualsUsed <- which(rownames(genonum) %in% as.character(envdata[,1]))
genome <- genonum[individualsUsed, ]

annotatedProbes <- read.table("annotatedExpressionData.txt", sep = "\t", header = TRUE)           # Read the expression data
expressiondata <- apply(as.matrix(annotatedProbes[,5:ncol(annotatedProbes)]), 2, log)
expressiondata <- normalizeBetweenArrays(expressiondata)

indIdx <- NULL
for(name in rownames(envdata)){ indIdx <- c(indIdx, grep(name, colnames(annotatedProbes))) }

colnames(annotatedProbes)[indIdx] <- as.character(envdata[,1])                                    # Transform the GSM into KT identifiers
genome <- genome[colnames(annotatedProbes)[indIdx],]                                              # Genome structure, now every row in genome has a matching column in annotatedProbes

mapQTL <- function(name = "AT3G16770_at", plot = TRUE){                                           # Map QTLs by gene name, the default should show a nice cis-eQTL
  idx <- which(annotatedProbes[,1] == name)
  lodMatrix <- NULL
  for(id in idx){
    phenotype <- 
    lod <- NULL
    for(x in 1:ncol(genome)){ lod <- c(lod, -log10(anova(lm(expressiondata[id,] ~ as.numeric(envdata[,2]) + genome[,x]))[[5]][2])) }
    if(plot) plot(lod, t='l')
    lodMatrix <- rbind(lodMatrix, lod)
  }
  rownames(lodMatrix) <- NULL
  colnames(lodMatrix) <- colnames(genome)
  invisible(list(annotation = annotatedProbes[idx,1:4], phenotypes = expressiondata[idx,], lodscores = lodMatrix))
}

lodscores <- mapQTL()

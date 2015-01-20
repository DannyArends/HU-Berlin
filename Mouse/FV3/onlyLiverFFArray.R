# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

library(preprocessCore)                                                                                               # preprocessCore - For quantile normalisation

setwd("E:/Mouse/RNA/FV3")

rawdata  <- read.table("RawData/expressions.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)              # Load the expression data
arrays   <- read.table("Annotation/arrays.txt", sep="\t", header=TRUE, row.names=1)                                   # Arrays annotation
probes   <- read.table("Annotation/illumina.txt", sep="\t", header=TRUE)                                              # Probe annotation

rawdata <- rawdata[, grep("1562232028", colnames(rawdata))]                           # Liver FF array
arrays <- arrays[grep("1562232028", rownames(arrays)),]                               # Liver FF array

# Preprocessing
rawdata[,rownames(arrays)] <- log2(rawdata[,rownames(arrays)])                              # Log2 transformation
rawdata[,rownames(arrays)] <- normalize.quantiles(as.matrix(rawdata[,rownames(arrays)]))    # Quantile normalisation

write.table(cbind(decoder_ID = rownames(rawdata), rawdata), "Analysis/NormDataLiverArray.txt", sep = "\t", quote=FALSE, row.names=FALSE)

boxplot(rawdata[,rownames(arrays)])

source("D:/Github/HU-Berlin/Mouse/FV3/createAnnotation.R")

if(!file.exists("Analysis/geneexpressionLiverArray.txt")){
  alldata <- NULL
  cnt <- 1
  ensgenes <- as.character(unique(annotationmatrix[,"ensembl_gene_id"]))
  for(x in ensgenes){
    annotsubset <- annotationmatrix[which(annotationmatrix[,"ensembl_gene_id"] == x),]
    cat(paste0(cnt, "/", length(ensgenes), ", Gene:"), x, ", Probes:", dim(annotsubset)[1], "\n")
    probeinformation <- NULL
    for(y in as.character(annotsubset[,"ProbeName"])){
      idx <- which(rownames(rawdata) == y)
      if(length(idx) > 0) probeinformation <- rbind(probeinformation, rawdata[idx,])
    }
    
    if(!is.null(probeinformation)){
      annotsubset <- annotsubset[annotsubset[,"ProbeName"]%in% rownames(probeinformation), ]
      alldata <- rbind(alldata, cbind(annotsubset, probeinformation))
    }
    cnt <- cnt + 1
  }
  write.table(alldata, file="Analysis/geneexpressionLiverArray.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading gene expression data from disk\n")
  alldata <- read.table("Analysis/geneexpressionLiverArray.txt", sep="\t", header=TRUE, check.names=FALSE)
}

B6 <- rownames(arrays[arrays[,"line"] == "B6",])
BFMI <- rownames(arrays[arrays[,"line"] == "BFMI860",])

pvalues <- apply(alldata[,c(B6,BFMI)], 1, function(x){
    t.test(x[1:3], x[4:6])$p.value
  }
)

alldata <- cbind(alldata, tTest = pvalues)

cat(as.character(alldata[which(alldata[,"tTest"] < 0.005), "ensembl_gene_id"]), sep = "\n")

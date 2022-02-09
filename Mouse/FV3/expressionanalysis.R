# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

library(preprocessCore)                                                                                               # preprocessCore - For quantile normalisation

setwd("D:/Edrive/Mouse/RNA/FV3")

rawdata  <- read.table("RawData/Old/expressions.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)              # Load the expression data
arrays   <- read.table("Annotation/arrays.txt", sep="\t", header=TRUE, row.names=1)                                   # Arrays annotation
probes   <- read.table("Annotation/illumina.txt", sep="\t", header=TRUE)                                              # Probe annotation


# QC of our probes
boxplot(rawdata[,rownames(arrays)])

badarrays <- c("1562232038_E", "1562232037_A","1562232040_C", "1740174012_D", "1754402037_A")

rawdata <- rawdata[, -which(colnames(rawdata) %in% badarrays)]                           # Misbehaving arrays
arrays <- arrays[-which(rownames(arrays) %in% badarrays),]                               # Misbehaving arrays

# Preprocessing
rawdata[,rownames(arrays)] <- log2(rawdata[,rownames(arrays)])                              # Log2 transformation
rawdata[,rownames(arrays)] <- normalize.quantiles(as.matrix(rawdata[,rownames(arrays)]))    # Quantile normalisation

write.table(cbind(decoder_ID = rownames(rawdata), rawdata), "Analysis/NormData.txt", sep = "\t", quote=FALSE, row.names=FALSE)

# QC of our probes
boxplot(rawdata[,rownames(arrays)])
corM <- cor(rawdata[,rownames(arrays)], method="spearman")
rownames(corM) <- paste0(arrays[,"tissue"],"_",arrays[,"diet"],"_",arrays[,"line"])
colnames(corM) <- arrays[,"tissue"]
heatmap(corM)

source("D:/Github/HU-Berlin/Mouse/FV3/createAnnotation.R")

if(!file.exists("Analysis/geneexpression.txt")){
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
  write.table(alldata, file="Analysis/geneexpression.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading gene expression data from disk\n")
  alldata <- read.table("Analysis/geneexpression.txt", sep="\t", header=TRUE)
}

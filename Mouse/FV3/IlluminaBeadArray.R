# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

source("http://www.bioconductor.org/biocLite.R")
biocLite(c("beadarray", "limma", "GEOquery", "illuminaMousev1.db", "illuminaMousev2.db", "BeadArrayUseCases"))
biocLite(c("GOstats", "GenomicRanges", "Biostrings"))

library(beadarray)

setwd("E:/Mouse/RNA/FV3")
arrays    <- read.table("Annotation/arrayannotation.txt", sep="\t", header=TRUE, row.names=1)                                  # Arrays annotation

# Extract data from bead arrays

for(beadarray in unique(arrays[,"array"])){
  metrics <- read.table(paste0("RawData/", beadarray, "/Metrics.txt"), sep = "\t", header = TRUE , as.is = TRUE)
  snr <- metrics$P95Grn / metrics$P05Grn
  labs <- paste(metrics[, 2], metrics[, 3], sep = "_")
  par(mai = c(1.5, 0.8, 0.3, 0.1))
  plot(1:12, snr , pch = 19, ylab = "P95/P05", xlab = "", main = "Signal-to-noise ratio for our data", axes = FALSE , frame.plot= TRUE)
  axis(2)
  axis(1, 1:12, labs , las = 2)
  
  mydata <- readIllumina(dir = paste0("RawData/", beadarray), useImages = FALSE, illuminaAnnotation = "Mouse")
  boxplot(mydata, transFun = logGreenChannelTransform ,col= "green", ylab = expression(log[2](intensity)), las = 2, outline = FALSE , main = paste0(sectionNames(mydata)[1], " MAQC data"))
  datasumm.c <- summarize(BLData = mydata, useSampleFac = TRUE, sampleFac=as.character(unlist(mydata@sectionData$SampleGroup)))
  write.table(exprs(datasumm.c), file = paste0("Intermediate/summary",beadarray,".txt"), sep="\t", row.names=TRUE)
}

# Reload the extracted files and combine them into a single expression matrix

alldata <- vector("list", length(unique(arrays[,"array"])))
cnt <- 1
for(beadarray in unique(arrays[,"array"])){
  datasumm <- read.table(paste0("summary",beadarray,".txt"), sep="\t", check.names=FALSE)
  alldata[[cnt]] <- datasumm
  cnt <- cnt+1
}
allprobes <- unique(unlist(lapply(alldata,function(x){rownames(x)})))

alldata <- matrix(NA, length(allprobes), length(rownames(arrays)), dimnames=list(allprobes, rownames(arrays)))
for(beadarray in unique(arrays[,"array"])){
  datasumm <- read.table(paste0("summary",beadarray,".txt"), sep="\t", check.names=FALSE)
  for(x in 1:ncol(datasumm)){
    alldata[rownames(datasumm), colnames(datasumm)[x]] <- datasumm[,x]
  }
}
write.table(alldata, file="Intermediate/summaryArrays.txt", sep ="\t", row.names=TRUE)

# After creating the files (also create the annotation matrix using createAnnotation.R)

library(preprocessCore)                                                                                               # preprocessCore - For quantile normalisation

setwd("E:/Mouse/RNA/FV3")
arrays <- read.table("Annotation/arrayannotation.txt", sep="\t", header=TRUE, row.names=1)                                  # Arrays annotation
annotationmatrix <- read.table("Annotation/probeannotation.txt", sep="\t", header=TRUE)                                     # Probe annotation
rawdata <- read.table("Intermediate/summaryArrays.txt", sep = "\t", header=TRUE, check.names=FALSE)                         # Raw data extracted from the arrays

badarrays <- c("1562232038_E", "1562232037_A","1562232040_C", "1740174012_D", "1754402037_A")

rawdata <- rawdata[, -which(colnames(rawdata) %in% badarrays)]                           # Remove the misbehaving arrays
arrays <- arrays[-which(rownames(arrays) %in% badarrays),]                               # Remove the misbehaving arrays

# Preprocessing
normdata <- normalize.quantiles(as.matrix(rawdata))                                                                         # Quantile normalisation
colnames(normdata) <- colnames(rawdata)
rownames(normdata) <- rownames(rawdata)

normdata <- normdata[which(rownames(normdata) %in% annotationmatrix[,"ProbeName"]),]

alldata <- cbind(annotationmatrix[match(rownames(normdata), annotationmatrix[,"ProbeName"]),], normdata)[1:10,]
write.table(alldata, file="Analysis/geneexpression.txt", sep="\t", row.names=FALSE)

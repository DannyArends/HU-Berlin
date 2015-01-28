source("http://www.bioconductor.org/biocLite.R")
biocLite(c("beadarray", "limma", "GEOquery", "illuminaMousev1.db", "illuminaMousev2.db", "BeadArrayUseCases"))
biocLite(c("GOstats", "GenomicRanges", "Biostrings"))

library(beadarray)

setwd("E:/Mouse/RNA/2009 FV3 Muscle Tissue/")

metrics <- read.table("5022073007/Metrics.txt", sep = "\t", header = TRUE , as.is = TRUE)
snr <- metrics$P95Grn / metrics$P05Grn
labs <- paste(metrics[, 2], metrics[, 3], sep = "_")
par(mai = c(1.5, 0.8, 0.3, 0.1))
plot(1:12, snr , pch = 19, ylab = "P95/P05", xlab = "", main = "Signal-to-noise ratio for our data", axes = FALSE , frame.plot= TRUE)
axis(2)
axis(1, 1:12, labs , las = 2)

### Read in the Illumina array data on the muscle tissue
mydata <- readIllumina(dir = "5022073007", useImages = FALSE, illuminaAnnotation = "Mousev1")

### Show some raw data
mydata@sectionData
sectionNames(mydata)
numBeads(mydata)
head(mydata[[1]])

### Show some raw bead level data
getBeadData(mydata, array = 1, what = "Grn")[1:5]
getBeadData(mydata, array = 1, what = "GrnX")[1:5]
getBeadData(mydata, array = 1, what = "ProbeID")[1:5]

### Boxplot the different arrays (log2 scale)
boxplot(mydata, transFun = logGreenChannelTransform ,col= "green", ylab = expression(log[2](intensity)), las = 2, outline = FALSE , main = "HT -12 MAQC data")

### Plot intensities per array
par(mfrow = c(6, 2))
par(mar = c(1, 1, 1, 1))
imageplot(mydata, array = 1, high = "darkgreen", low = "lightgreen",zlim =c(4, 10), main = sectionNames(mydata)[1])

### Plot outliers per array
outlierplot(mydata, array = 1, main = paste(sectionNames(mydata)[1], "outliers"))

### Remove spatial artefacts
#for(x in 1:nrow(mydata@sectionData$Targets)){
#  BASHoutput <- BASH(mydata, array = x)
#  mydata <- setWeights(mydata, wts = BASHoutput$wts, array = x)
#}

datasumm <- summarize(BLData = mydata, useSampleFac = TRUE)

grnchannel <- new("illuminaChannel", transFun = greenChannelTransform  , outlierFun = illuminaOutlierMethod , exprFun = function(x) mean(x, na.rm = TRUE), varFun = function(x) sd(x,na.rm = TRUE), channelName = "G")
grnchannel.log <- new("illuminaChannel", transFun = logGreenChannelTransform , outlierFun = illuminaOutlierMethod , exprFun = function(x) mean(x, na.rm = TRUE), varFun = function(x) sd(x,na.rm = TRUE), channelName = "G")

datasumm.all <- summarize(BLData = mydata, useSampleFac = FALSE , channelList = list(grnchannel, grnchannel.log))


exprs(datasumm)[1:10,]
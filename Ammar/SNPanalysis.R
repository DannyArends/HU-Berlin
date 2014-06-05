# SNPanalysis.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2014
# first written May, 2014
# 
# Script for SNP analysis in Sudanese dairy cows

setwd("d:/ammar")
library(randomForest)
rawdata <- read.table("genotypes.txt")

numdata <- NULL
for(x in 1:ncol(rawdata)){
  numdata <- cbind(numdata, as.numeric(rawdata[,x]))
}

rownames(numdata) <- rownames(rawdata)
colnames(numdata) <- colnames(rawdata)

groups <- c("WN", "KN", "BT", "IRS", "ALG")
colz <- c("green", "orange", "purple", "blue", "yellow")

clusters <- hclust(dist(numdata))
orderinheatmap <- rownames(rawdata)[clusters$order]

rc <- rep("black", nrow(numdata))
ccol <- 1
for(x in groups){
  ids <- grep(x, orderinheatmap)
  rc[ids] <- colz[ccol]
  ccol <- ccol + 1
}

op <- par(cex.axis = 1.5)
op <- par(cex = 1.5)
png(file="heatmap.png", width=3000, height=3000)
dendro <- heatmap(numdata, Colv=NA, RowSideColors = rc, keep.dendro=TRUE)
dev.off()

op <- par(cex = 0.4)
png(file="dendrogram.png", width=2800, height=1000)
plot(dendro$Rowv)
dev.off()

# randomForest approach to determine selection criteria
groupVector <- rep(0, nrow(numdata))
for(x in groups){ groupVector[grep(x, rownames(numdata))] <- x }

ammardata <- cbind(numdata, groupVector)
colnames(ammardata) <- c(colnames(numdata), "Race")
ammar.rf <- randomForest(Race ~ ., data=ammardata, importance=TRUE, proximity=TRUE)
ammar.rf

round(importance(ammar.rf), 1)                                                                          ## Variable importance:
ammar.mds <- cmdscale(ammar.rf$proximity, eig=TRUE)                                                     ## Do MDS on 1 - proximity
colorz <- c("green", "orange", "purple", "blue", "yellow")[as.numeric(as.factor(ammardata[,"Race"]))]   ## Colors for the plot

png(file="randomForest.png", width=2000, height=2000)
op <- par(cex = 2.5)
plot(cbind(ammar.mds$points), cex=2.0, pch=19, col=colorz, main="Ammar Data: MDS of Proximity (RandomForest)", xlab="Var1", ylab="Var2")
dev.off()
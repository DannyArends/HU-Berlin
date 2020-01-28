# Source some functions
source("D:/Ddrive/Github/HU-Berlin/Goat/Africa/functions.R")

setwd("D:/Edrive/Goat/")
data.all <- read.table(file = "Africa.Swiss.Adaptmap.merged.matrix.noDuplicates.txt", sep = "\t", colClasses = "character")
data.num <- read.table(file = "Africa.Swiss.Adaptmap.merged.num.matrix.noDuplicates.txt", sep = "\t", colClasses = "character")
snp.annot <- read.table(file = "Africa.Swiss.Adaptmap.snp.annotation.txt", sep = "\t", colClasses = "character")
sample.annot <- read.table(file = "Africa.Swiss.Adaptmap.samples.annotation.noDuplicates.txt", sep = "\t", quote="\"", colClasses = "character")
dist.mat <- read.table("Africa.Swiss.Adaptmap.distancematrix.noDuplicates.txt", sep = "\t")

# Fix annotation errors
sample.annot[which(sample.annot[, "Breed"] == "Small East Africa"), "Breed"] <- "Small East African"
# Cluster to find the outliers
clusters <- hclust(as.dist(dist.mat))

# remove outliers and missing breeds, and cluster again
outliers <- which(cutree(clusters,2) == 2)
missingBreed <- which(is.na(sample.annot[, "Breed"]))
mixes <- grep(" x ", sample.annot[, "Breed"])
mixes <- c(mixes, grep("Admixed", sample.annot[, "Breed"]))
related <- tooCloselyRelated(dist.mat)

toremove <- c(outliers, missingBreed, mixes, related)

clusters <- hclust(as.dist(dist.mat[-toremove, -toremove]))

dendrogram <- as.dendrogram(clusters)   

continents <- 1:6
names(continents) <- c("Africa", "Europe", "N America", "Oceania", "S America", "W Asia")

svg("Africa.Swiss.Adaptmap.dendrogram.identifiers.svg", width = 300, height = 7, pointsize = 6)
  dendroColor <- dendrapply(dendrogram, labelCol)
  plot(dendroColor, cex=0.7)
dev.off()

svg("Africa.Swiss.Adaptmap.dendrogram.breeds.svg", width = 300, height = 7, pointsize = 6)
  dendroColor <- dendrapply(dendrogram, labelCol, LblAsBreed = TRUE)
  plot(dendroColor, cex=0.7)
dev.off()

maxDistances <- apply(dist.mat,1,max)
hist(maxDistances, breaks = 100)

#
# Continue with only the full data by removing animals as above
#

data.num <- data.num[,-toremove]
data.all <- data.all[,-toremove]
sample.annot <- sample.annot[-toremove,]

marker.missing <- apply(apply(data.num,1,is.na),2,sum)
sample.missing <- apply(apply(data.num,2,is.na),2,sum)

m10pct <- which(marker.missing > 0.1 * ncol(data.num))
s10pct <- which(sample.missing > 0.1 * nrow(data.num))

data.num <- data.num[-m10pct, -s10pct]
data.all <- data.all[-m10pct, -s10pct]
sample.annot <- sample.annot[-s10pct,]
snp.annot <- snp.annot[-m10pct,]

write.table(data.all, file = "merged.filtered.matrix.txt", sep = "\t", quote=FALSE)
write.table(data.num, file = "merged.filtered.num.matrix.txt", sep = "\t", quote=FALSE)
write.table(sample.annot, file = "merged.samples.annotation.txt", sep = "\t", quote=FALSE)
write.table(snp.annot, file = "merged.snp.annotation.txt", sep = "\t", quote=FALSE)

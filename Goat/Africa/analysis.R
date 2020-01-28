# Source some functions
source("D:/Ddrive/Github/HU-Berlin/Goat/Africa/functions.R")

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
map <- read.table("FilteredLocationLWLT01.txt", sep = "\t", header=TRUE, row.names=1)
map <- chromosomeOrder(map)

setwd("D:/Edrive/Goat/")

data.all <- read.table(file = "merged.filtered.matrix.txt", sep = "\t", colClasses = "character")
data.num <- read.table(file = "merged.filtered.num.matrix.txt", sep = "\t", colClasses = "character")
snp.annot <- read.table(file = "merged.snp.annotation.txt", sep = "\t", colClasses = "character")
sample.annot <- read.table(file = "merged.samples.annotation.txt", sep = "\t", quote="\"", colClasses = "character")

data.num <- data.num[which(rownames(data.num) %in% rownames(map)),]
map <- map[which(rownames(map) %in% rownames(data.num)),]

data.part <- apply(data.num[rownames(map),], 2, as.numeric)
marker.cor <- cor(t(data.part), use = "pair")
write.table(marker.cor, file = "merged.filtered.marker.correlation.matrix.txt", sep = "\t", quote=FALSE)

marker.cor.round <- apply(marker.cor, 2, round, 2)
write.table(marker.cor.round, file = "merged.filtered.marker.correlation.2dig.matrix.txt", sep = "\t", quote=FALSE)

# Per breed analysis to figure out 'fixed' SNP alleles based on our breed data
all.breeds <- unique(sample.annot[,"Breed"])
good.breeds <- names(which(table(sample.annot[,"Breed"]) >= 40))

sample.annot <- sample.annot[which(sample.annot[,"Breed"] %in% good.breeds),]
data.num <- data.num[,rownames(sample.annot)]

consistentGTs <- matrix(NA, nrow(data.num), length(good.breeds), dimnames = list(rownames(data.num), good.breeds))
for(breed in good.breeds){
  ind.names <- rownames(sample.annot)[sample.annot[, "Breed"] == breed]
  tbl <- apply(data.num[,ind.names],1,table)
  idx <- which(unlist(lapply(tbl, length)) == 1)
  consistentGTs[idx, breed] <- unlist(lapply(tbl[idx], names))
  cat("Done", breed, "\n")
}

consistentGTs <- consistentGTs[which(apply(consistentGTs, 1, function(x){sum(is.na(x)) != length(x)})),]


plot(as.numeric(table(sample.annot[,"Breed"])), apply(consistentGTs, 2, function(x){sum(!is.na(x))}))
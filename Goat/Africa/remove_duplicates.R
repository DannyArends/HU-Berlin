setwd("D:/Edrive/Goat/")

data.all <- read.table(file = "Africa.Swiss.Adaptmap.merged.matrix.txt", sep = "\t", colClasses = "character")
data.num <- read.table(file = "Africa.Swiss.Adaptmap.merged.num.matrix.txt", sep = "\t", colClasses = "character")
sample.annot <- read.table(file = "Africa.Swiss.Adaptmap.breed.annotation.txt", sep = "\t", quote="\"")
dist.matFull <- read.table("Africa.Swiss.Adaptmap.distancematrix.txt", sep = "\t")

# Read in the distance matrix and remove the lower triangle and diagonal
dist.mat <- read.table("Africa.Swiss.Adaptmap.distancematrix.txt", sep = "\t")
dist.mat[lower.tri(dist.mat, diag = TRUE)] <- NA

duplicates <- c()
for(x in 1:nrow(dist.mat)){
  idx <- which(dist.mat[x, ] == 0)
  if (length(idx) > 0) {
    for (i in idx) {
      n1 <- rownames(dist.mat)[x]
      n2 <- rownames(dist.mat)[i]
      cat("Duplicate sample:", n1, "with", n2,"\n")
      duplicates <- rbind(duplicates, c(n1,n2))
    }
  }
}

dist.idx <- which(rownames(dist.matFull) %in% duplicates[,2])
dist.matFull <- dist.matFull[-dist.idx, -dist.idx]
write.table(dist.matFull, file = "Africa.Swiss.Adaptmap.distancematrix.noDuplicates.txt", sep = "\t", quote=FALSE)

sample.annot <- sample.annot[-which(rownames(sample.annot) %in% duplicates[,2]),]
write.table(sample.annot, file = "Africa.Swiss.Adaptmap.samples.annotation.noDuplicates.txt", sep = "\t", quote=FALSE)

data.num <- data.num[, -which(colnames(data.num) %in% duplicates[,2])]
write.table(data.num, file = "Africa.Swiss.Adaptmap.merged.num.matrix.noDuplicates.txt", sep = "\t", quote=FALSE)

data.all <- data.all[, -which(colnames(data.all) %in% duplicates[,2])]
write.table(data.all, file = "Africa.Swiss.Adaptmap.merged.matrix.noDuplicates.txt", sep = "\t", quote=FALSE)

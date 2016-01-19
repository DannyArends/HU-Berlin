# Analysis of horse SNP chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2016
# first written Jan, 2016

toNumeric <- function(genotypes){
  numericG <- apply(genotypes, 1, function(x){
    geno <- table(unlist(strsplit(as.character(x),"")))
    a1 <- paste0(names(geno)[1],names(geno)[1])
    a2 <- paste0(sort(c(names(geno)[1],names(geno)[2])),collapse="")
    a3 <- paste0(names(geno)[2],names(geno)[2])
    ngeno <- rep(NA,length(x))
    ngeno[x == a1] <- 1
    ngeno[x == a2] <- 2
    ngeno[x == a3] <- 3
    return(ngeno)
  })
  rownames(numericG) <- colnames(genotypes)
  return(t(numericG))
}

setwd("E:/Horse/DNA/")
if(!file.exists("Equine60k/input/cleaned_arabian_petersen.txt")){
  # Create a matrix of genotypes
  raw_data <- read.csv("Petersen2013/ArabianHorses.csv", header = FALSE, sep = ",")
  samples <- as.character(unique(raw_data[,2]))
  markers <- unique(raw_data[,1])

  mdata <- matrix(NA, length(markers), length(samples), dimnames=list(markers, samples))
  for(x in 1:nrow(raw_data)){
    if(!is.nan(raw_data[x, 7]) && raw_data[x, 7] >= 0.15) mdata[as.character(raw_data[x, 1]), as.character(raw_data[x, 2])] <- paste0(raw_data[x, 5], raw_data[x, 6])
  }

  rownames(mdata) <- gsub("-", "_", rownames(mdata))
  write.table(mdata, file="Equine60k/input/cleaned_arabian_petersen.txt", sep = "\t")
}else{
  mdata <- read.table("Equine60k/input/cleaned_arabian_petersen.txt", sep = "\t")
}

# Read our data
map <- read.table(file="Equine60k/input/cleaned_map.txt", sep = "\t")                                      # Save the clean map to disk
genotypes <- read.table(file="Equine60k/input/cleaned_genotypes.txt", sep = "\t")                          # Save the clean genotypes to disk


mdata <- mdata[which(rownames(mdata) %in% rownames(genotypes)),]  # 28K markers are shared
genotypes <- genotypes[which(rownames(genotypes) %in% rownames(mdata)),]  # 28K markers are shared


combined <- cbind(genotypes, mdata)

numGeno <- t(toNumeric(combined))


## Some basic plots of all the individuals relatedness
dendrogram <- as.dendrogram(hclust(dist(numGeno, method = "manhattan")))         #TODO: perhaps add the reference horse

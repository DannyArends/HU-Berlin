# Analysis of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016

setwd("E:/Goat/DNA/SihamAnalysis")

snpdata <- read.table("filtered_snps.txt", sep="\t", check.names=FALSE, colClasses="character")
snpinfo <- read.table("snpinfo.txt", sep="\t")
samples <- read.table("sampleinfo.txt", sep="\t")

snpAlleles <- lapply(strsplit(as.character(snpinfo[,"allele"]), ""), "[", c(1,3))

if(!file.exists("filtered_snps_numeric.txt")){
  numsnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames = list(rownames(snpdata), colnames(snpdata)))
  for(x in 1:length(snpAlleles)) {
    if(!is.na(snpinfo[x, "reference"]) && snpAlleles[[x]][1] !=  snpinfo[x, "reference"]){  # C/T while reference is T, so flip it around
      snpAlleles[[x]] <- snpAlleles[[x]][2:1]
    }

    g1 <- paste(snpAlleles[[x]][1], snpAlleles[[x]][1],sep="")
    g2a <- paste(snpAlleles[[x]][1], snpAlleles[[x]][2],sep="")
    g2b <- paste(snpAlleles[[x]][2], snpAlleles[[x]][1],sep="")
    g3 <- paste(snpAlleles[[x]][2], snpAlleles[[x]][2],sep="")
    if(!all(snpdata[x,] %in% c(g1,g2a,g2b,g3, NA))) stop("Nope")
    numsnpdata[x, which(snpdata[x, ] == g1)] <- 1
    numsnpdata[x, which(snpdata[x, ] == g2a)] <- 2
    numsnpdata[x, which(snpdata[x, ] == g2b)] <- 2
    numsnpdata[x, which(snpdata[x, ] == g3)] <- 3
  }

  write.table(numsnpdata, "filtered_snps_numeric.txt", sep="\t", quote=FALSE)
}else{
  numsnpdata <- read.csv("filtered_snps_numeric.txt", sep="\t", check.names=FALSE)
}

# Add the reference animal to the dataset, since reference is always coded as 1
#numsnpdata <- cbind(numsnpdata, Reference = 1)
clustering <- hclust(dist(t(numsnpdata), method = "manhattan"))

plot(as.phylo(clustering), type="r")  # Or a rooted dendrogram plot: plot(root(as.phylo(clustering),"Reference"), type="r")
plot(clustering, hang = -1)
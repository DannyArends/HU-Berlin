library(heterozygous)
library(genetics)

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")

#Data is filtered for: GenTrain > 0.6, MAF > 0.05, unknown positions, 5% call rate 
snpdata    <- read.table("filtered_snps_NO_DN2.txt", sep="\t", check.names=FALSE, colClasses="character")
snpinfo    <- read.table("snpinfo.txt", sep="\t", na.strings=c("", "NA", "N.D."))
samples    <- read.table("sampleinfo.txt", sep="\t")
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`

dim(snpdata)

# Filter MAF more aggresively <= 0.1
allelefreq <- apply(snpdata, 1 , function(x){
  tbl <- table(unlist(lapply(x, strsplit, "")))
  min(tbl / sum(tbl))
})

snpinfo <- cbind(snpinfo, MAF=allelefreq)

# Prune MAP for markers in LD R2 > 0.5, keep the ones with the highest MAF
map <- snpinfo[,c("Chr", "Position", "MAF")]
map.pruned <- map

genotypes <- snpdata
genotypes.pruned <- genotypes
markers <- rownames(genotypes.pruned)

x <- 1
while(x < length(markers)) {
  mName <- markers[x]
  mChr <- as.character(map[mName, "Chr"])
  mPos <- as.numeric(map[mName, "Position"])

  nearby <- which(as.character(map.pruned[,"Chr"]) == mChr & as.numeric(map.pruned[, "Position"]) > (mPos - 1000000) &  as.numeric(map.pruned[, "Position"]) < (mPos + 1000000))
  nearby <- rownames(map.pruned)[nearby]

  mGeno <- genotype(as.character(genotypes.pruned[mName,]), sep = "")
  locOfMarker <- which(nearby == mName)
  if(locOfMarker > 1) nearby <- nearby[-(1:locOfMarker-1)]

  LDs <- rep(NA, length(nearby))
  names(LDs) <- nearby
  for(y in nearby){
     LDs[y] <- round(LD(mGeno, genotype(as.character(genotypes.pruned[y,]), sep = ""))$"R^2",2)
  }
  inLD <- names(which(LDs > 0.5))
  if(length(inLD) > 1){
    MAFs <- map.pruned[inLD, ]
    bestSNP <- rownames(MAFs)[which.max(MAFs[, "MAF"])]
    inLD <- inLD[-which(inLD == bestSNP)]
    toPrune <- which(rownames(genotypes.pruned) %in% inLD)
    if(length(toPrune) > 0){
      genotypes.pruned <- genotypes.pruned[-toPrune, ]
      map.pruned <- map.pruned[-toPrune, ]
      markers <- rownames(genotypes.pruned)
    }
  }
  cat("Done", x, "/", nrow(genotypes), "/", nrow(genotypes.pruned), "==", nrow(map.pruned), "\n")
  x <- (x + 1)
}

write.table(genotypes.pruned, file="LDpruned_genotypes.txt", sep = "\t")
write.table(map.pruned, file="LDpruned_map.txt", sep = "\t")

snpinfo <- snpinfo[rownames(map.pruned),]

snpAlleles <- lapply(strsplit(as.character(snpinfo[,"allele"]), ""), "[", c(1,3))

numsnpdata <- matrix(NA, nrow(genotypes.pruned), ncol(genotypes.pruned), dimnames = list(rownames(genotypes.pruned), colnames(genotypes.pruned)))
for(x in 1:length(snpAlleles)) {
  if(!is.na(snpinfo[x, "reference"]) && snpAlleles[[x]][1] !=  snpinfo[x, "reference"]){  # C/T while reference is T, so flip it around
    snpAlleles[[x]] <- snpAlleles[[x]][2:1]
  }

  g1 <- paste(snpAlleles[[x]][1], snpAlleles[[x]][1],sep="")
  g2a <- paste(snpAlleles[[x]][1], snpAlleles[[x]][2],sep="")
  g2b <- paste(snpAlleles[[x]][2], snpAlleles[[x]][1],sep="")
  g3 <- paste(snpAlleles[[x]][2], snpAlleles[[x]][2],sep="")
  if(!all(genotypes.pruned[x,] %in% c(g1,g2a,g2b,g3, NA))) stop("Nope")
  numsnpdata[x, which(genotypes.pruned[x, ] == g1)] <- 1
  numsnpdata[x, which(genotypes.pruned[x, ] == g2a)] <- 2
  numsnpdata[x, which(genotypes.pruned[x, ] == g2b)] <- 2
  numsnpdata[x, which(genotypes.pruned[x, ] == g3)] <- 3
}

write.table(numsnpdata, file="LDpruned_num_genotypes.txt", sep = "\t")


setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
numsnpdata <- read.table(file="LDpruned_num_genotypes.txt", sep = "\t")
map.pruned <- read.table(file="LDpruned_map.txt", sep = "\t")
samples    <- read.table("sampleinfo.txt", sep="\t")
snpinfo    <- read.table("snpinfo.txt", sep="\t", na.strings=c("", "NA", "N.D."))


### Principal component analysis
misdata <- which(apply(numsnpdata, 1, function(x){any(is.na(x))}))
numsnppca <- numsnpdata[-misdata,]

pcares <- prcomp(t(numsnppca), scale=TRUE)
groups <- samples[colnames(numsnppca), "Breed"]

sumpca <- summary(pcares)

pIC <- pcares$x[,1:10]
npIC <- apply(pIC, 2, function(x){
  r <- (max(x) - min(x))
  (((x - min(x)) / r) - 0.5)
})


pca1 <- paste0("(", round(sumpca$importance[2,1] * 100,1), "%", " variance explained)")
pca2 <- paste0("(", round(sumpca$importance[2,2] * 100,1), "%", " variance explained)")
pca3 <- paste0("(", round(sumpca$importance[2,3] * 100,1), "%", " variance explained)")



 # Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}

# Variable correlation/coordinates
var.coord <- t(apply(pcares$rotation, 1, var_cor_func, pcares$sdev))
head(var.coord[, 1:4])
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))


### Look into the rotations
importantSNPs <-  names(which(var.contrib[,1] >= 0.02))
lessimpSNPs <-  names(which(var.contrib[,1] >= 0.02))

#snpinfo <- snpinfo[-which(!rownames(snpinfo) %in% rownames(pcares$rotation)), ]

SNPsPCA1 <- snpinfo[importantSNPs, c("Chr", "Position")]
SNPlPCA1 <- snpinfo[lessimpSNPs, c("Chr", "Position")]



chromosomes <- c(as.character(1:29),"X")
ymax <- max(snpinfo[,"Position"])

op <- par(mar=c(5, 4, 1, 2) + 0.1)

#png("PCAplotLocations.png", width=1200, height=600, res=300, pointsize = 5)
plot(x=c(0,length(chromosomes)), y = c(0, ymax), t='n', xaxt='n', yaxt='n', ylab="", xlab="Chromosome")
chrid <- 1
for(chr in chromosomes){
  allG <- snpinfo[snpinfo[,"Chr"] == chr, "Position"]
  chrM <- max(snpinfo[snpinfo[,"Chr"] == chr, "Position"])
  intG <- SNPsPCA1[which(SNPsPCA1[,"Chr"] == chr),"Position"]
  intP <- SNPlPCA1[which(SNPlPCA1[,"Chr"] == chr),"Position"]
  lines(x=c(chrid,chrid), y = c(0, chrM *0.9))
  points(x = rep(chrid, length(allG)), allG, pch="-", col=rgb(0.9, 0.9, 0.9), cex=2)
  points(x = rep(chrid, length(intP)), intP, pch="-", col=rgb(0.5, 0.5, 0.5), cex=2)
  points(x = rep(chrid, length(intG)), intG, pch="-", col="black", cex=2)
  chrid <- chrid + 1
}
axis(1, at = 1:length(chromosomes), chromosomes, lwd=0, lwd.tick=0.4)
axis(2, at = seq(0, ymax, 10000000), paste(seq(0, ymax, 10000000) / 1000000, "Mb"),las=2, lwd=0, lwd.tick=0.4)

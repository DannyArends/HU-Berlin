
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/incompatible.R")

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")

# Load the genotype call data
phased.vcf <- read.table(gzfile(paste0("TRDredo/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("TRDredo/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header

samples <- colnames(phased.vcf)[-c(1:9)]

phased.vcf[, samples] <- apply(phased.vcf[, samples],2,function(x){unlist(lapply(strsplit(x,":"),"[",1))})                   # Only the genotypes
rownames(phased.vcf) <- phased.vcf[,"ID"]
phased.vcf <- phased.vcf[,-which(colnames(phased.vcf) %in% c("ID", "QUAL","FILTER","INFO","FORMAT"))]                        # Remove unneeded columns
phased.vcf[1:10,1:10]

# Load in the phenotypes and create the pedigree file
phenotypes <- read.table("Phenotypes/allPhenotypes.txt",sep="\t",header=TRUE,na.strings=c(0,"-","NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.vcf)),]                            # We do not have genotypes for all individuals

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals

# Change the genotype coding
phased.num <- phased.vcf                                                                         # Copy
phased.num[, samples] <- apply(phased.num[, samples], 2, fromVCFcode.Num)                        # Change coding to A H B

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
# Load the Allele Transmission Biased regions
regionsMat <- read.table("TRDredo/mat0.01annot.txt", sep="\t", header=FALSE)
regionsPat <- read.table("TRDredo/pat0.01annot.txt", sep="\t", header=FALSE)

inMat <- c()
for(x in 1:nrow(regionsMat)){
  inside <- rownames(phased.num)[which(phased.num[, "CHROM"] == regionsMat[x, 3] & phased.num[, "POS"] >= regionsMat[x, 4] & phased.num[, "POS"] <= regionsMat[x, 6])]
  inMat <- c(inMat, inside)
}

inPat <- c()
for(x in 1:nrow(regionsPat)){
  inside <- rownames(phased.num)[which(phased.num[, "CHROM"] == regionsPat[x, 3] & phased.num[, "POS"] >= regionsPat[x, 4] & phased.num[, "POS"] <= regionsPat[x, 6])]
  inPat <- c(inPat, inside)
}

outside.num <- phased.num[which(!(rownames(phased.num) %in% inMat) & !(rownames(phased.num) %in% inPat)),]
inside.num <- phased.num[which(!rownames(phased.num) %in% rownames(outside.num)),]

doAF <- function(dataset, ind){
  rr <- apply(dataset[,F2], 1, function(x){
    aa <- length(which(x == -1)) * 2 + length(which(x == 0))
    bb <- length(which(x == 1)) * 2 + length(which(x == 0))
    c(aa,bb)
  })
  af <- apply(rr,2,function(x){ min(x) / sum(x) })
  return(af[af > 0.1])
}

afOut <- doAF(outside.num, F2)
afIn <- doAF(inside.num, F2)
op <- par(mfrow = c(2,1))
hist(afOut, main="MAF markers outside TRD regions", breaks = seq(0.1,0.5,0.025), xlab="Minor Allele Frequency")
hist(afIn, main="MAF markers inside TRD regions", breaks = seq(0.1,0.5,0.025), xlab="Minor Allele Frequency")

gts <- t(phased.vcf[, rownames(phenotypes)])
gtsinfo <- cbind(phenotypes[, c("Gen.", "sex", "Vater", "Mutter")], gts)
colnames(gtsinfo)[1:4] <- c("Generation", "Sex", "Father", "Mother")

map <- read.table("TRDredo/markerAnnotation.txt",sep="\t")
mmap <- t(map[colnames(gtsinfo),1:2])
colnames(mmap)[1:4] <- c("Generation", "Sex", "Father", "Mother")

gtsinfo <- rbind(mmap, gtsinfo)

write.table(t(gtsinfo), file="TRDredo/Phased_Genotypes_Info.txt", sep = "\t", quote=FALSE)

locusxdna <- read.table("TRDredo/genotypes.txt", sep="\t", check.names=FALSE)

gtsinfo <- rbind(gtsinfo[1:2,], B6N=NA, "BFMI860-12 (V2)"=NA, gtsinfo[3:nrow(gtsinfo),])

gtsinfo[3:nrow(gtsinfo), 5:ncol(gtsinfo)] <- t(locusxdna[colnames(gtsinfo[-c(1:2), -c(1:4)]), c(rownames(gtsinfo[-c(1:2), -c(1:4)]))])
gtsinfo <- t(gtsinfo)

order <- sort(gtsinfo[1,], index=TRUE,na.last=FALSE,method = "radix")$ix

write.table(gtsinfo[, order], file="TRDredo/Raw_Genotypes_Info.txt", sep = "\t", quote=FALSE,na="-")
mmap <- t(mmap)

# Load the Allele Transmission Biased regions
regionsMat <- read.table("TRDredo/mat0.01annot.txt", sep="\t", header=FALSE)
regionsPat <- read.table("TRDredo/pat0.01annot.txt", sep="\t", header=FALSE)


size <- 5000000

plot(c(1, max(map[,"Pos"],na.rm=TRUE)), c(0,20), main="Marker density", ylab="Chromosome", xlab="Postion (Mb)", t ='n', yaxt='n', xaxt='n')
cnt <- 1
nMmax <- c(0)
for(x in c(1:19,"X")){
  chrmap <- mmap[which(mmap[, "Chr"] == x),]
  chrmax <- max(map[which(map[, "Chr"] == x), "Pos"],na.rm=TRUE)#max(as.numeric(chrmap[, "Pos"]))
  for(s in seq(0, chrmax - (size/5), size / 5)){
    e <- s + (size/5)
    nM <- length(which(as.numeric(chrmap[, "Pos"]) > s & as.numeric(chrmap[, "Pos"]) < e))
    rect(s, cnt, e, cnt + nM / 100, col="gray")
    #text(s+ size/2, cnt, nM, cex=(69+nM) / 205)
  }
  cnt <- cnt+1
}

cnt <- 1
nMmax <- c(0)
for(x in c(1:19,"X")){
  chrmap <- mmap[which(mmap[, "Chr"] == x),]
  chrmax <- max(map[which(map[, "Chr"] == x), "Pos"],na.rm=TRUE)#max(as.numeric(chrmap[, "Pos"]))
  for(s in seq(0, chrmax - size, size)){
    e <- s + size
    nM <- length(which(as.numeric(chrmap[, "Pos"]) > s & as.numeric(chrmap[, "Pos"]) < e))
    if(nM > nMmax) nMmax <- nM
    text(s+ size/2, cnt-0.2, nM, cex=0.5, col="black")
  }
  cnt <- cnt+1
}

for(x in 1:nrow(regionsMat)){
  cnt <- regionsMat[x,3]
  s <- regionsMat[x,4]
  e <- regionsMat[x,6]
  cat(cnt, s, e, "\n")
  rect(s, cnt-0.05, e, cnt - 0.1, col="green", border = NA)
}

for(x in 1:nrow(regionsPat)){
  cnt <- regionsPat[x,3]
  s <- regionsPat[x,4]
  e <- regionsPat[x,6]
  cat(cnt, s, e, "\n")
  rect(s, cnt-0.05, e, cnt - 0.1, col="green", border = NA)
}

axis(2, at = 1:20, c(1:19,"X"), las=2)
axis(1, at = seq(0, 200000000, 10000000), seq(0, 200000000, 10000000) / 1000000, las=2)

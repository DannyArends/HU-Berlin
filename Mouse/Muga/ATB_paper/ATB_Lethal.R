#
# Analyse ATB in megaMuga data from multiple generations after beagle phasing
#

# Load the genotype call data

source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

setwd("E:/Mouse/DNA/MegaMuga/")

# Load the genotype call data
phased.vcf <- read.table(gzfile(paste0("Analysis/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("Analysis/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header

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
phased.AHB <- phased.vcf                                                                         # Copy
phased.AHB[, samples] <- apply(phased.AHB[, samples], 2, fromVCFcode.AHB)                        # Change coding to A H B

# QTL Analysis of generation 28
# Filter the genotypes for markers that are seggregating
phased.AHB <- phased.AHB[-which(lapply(apply(phased.AHB[, F2], 1, table),length) == 1),]

genotypes <- NULL
for(x in c(seq(1,19),"X","Y","M")) {  # Reorder the chromosomes and name it genotypes
  genotypes <- rbind(genotypes, phased.AHB[which(phased.AHB[,"CHROM"] == x), ])
}



# Load the Allele Transmission Biased regions
regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)

sift <- read.table("siftscores.txt", sep="\t")
sift <- sift[which(sift[,16] == "Deleterious" | sift[,20] == "Damaging"),]
sift <- sift[which(!duplicated(sift[,2])),]

nBFMI <- 0;sumBFMI <- 0;
nB6N <- 0;sumB6N <- 0;
for(x in 1:nrow(regionsMat)){
  chr <- regionsMat[x, "Chr"]
  sta <- regionsMat[x, "Start"]
  sto <- regionsMat[x, "Stop"]
  allele <- as.character(regionsMat[x, "Prefered.Allele"])
  snpL <- as.numeric(sift[,4])
  nSiftAlleles <- length(which(sift[,3] == chr & snpL > sta & snpL < sto))
  if(allele == "BFMI"){
    nBFMI <- nBFMI + 1
    sumBFMI <- sumBFMI + nSiftAlleles
  }
  if(allele == "B6N"){
    nB6N <- nB6N + 1
    sumB6N <- sumB6N + nSiftAlleles
  }

  cat(x, allele, nSiftAlleles, "\n")
}

cat("nBFMI:", nBFMI, "SUM:", sumBFMI, "\nnB6N", nB6N, "SUM:", sumB6N, "\n")


nBFMI <- 0;sumBFMI <- 0;
nB6N <- 0;sumB6N <- 0;
for(x in 1:nrow(regionsPat)){
  chr <- regionsPat[x, "Chr"]
  sta <- regionsPat[x, "Start"]
  sto <- regionsPat[x, "Stop"]
  allele <- as.character(regionsPat[x, "Prefered.Allele"])
  snpL <- as.numeric(sift[,4])
  nSiftAlleles <- length(which(sift[,3] == chr & snpL > sta & snpL < sto))
  if(allele == "BFMI"){
    nBFMI <- nBFMI + 1
    sumBFMI <- sumBFMI + nSiftAlleles
  }
  if(allele == "B6N"){
    nB6N <- nB6N + 1
    sumB6N <- sumB6N + nSiftAlleles
  }

  cat(x, allele, nSiftAlleles, "\n")
}

cat("nBFMI:", nBFMI, "SUM:", sumBFMI, "\nnB6N", nB6N, "SUM:", sumB6N, "\n")


map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
aa <- read.table("analysis/QTLresults_newModelQTL.txt",sep="\t")


xmap <- map[which(map[,1] == "X"),]

markers <- rownames(xmap[which(xmap[,2] > 50000000 & xmap[,2] < 80000000),])
markers <- markers[which(markers %in% rownames(aa))]

plot(aa[markers,"d70"])


allRegions <- rbind(cbind(regionsMat, origin = "M"), cbind(regionsPat, origin = "P"))
allRegions <- allRegions[with(allRegions, order(Chr, Start)), ]
rownames(allRegions) <- 1:nrow(allRegions)

## new idea compare 2 regions against each other:
results <- matrix(0, nrow(allRegions), nrow(allRegions))
results2 <- matrix(0, nrow(allRegions), nrow(allRegions))

sG <- c("AA", "AB", "BA", "BB")

for(y in 1:nrow(allRegions)){
  mR1 <- genotypes[as.character(allRegions[y,"Flanking.Marker"]), F2]
  chrR1 <- allRegions[y,"Chr"]
  for(x in 1:nrow(allRegions)) {
    mR2 <- genotypes[as.character(allRegions[x,"Flanking.Marker"]), F2]
    chrR2 <- allRegions[x,"Chr"]
    if(chrR1 == chrR2){
      results[y,x] <- NA
    }else{
      observed <- rep(0, 9)
      names(observed) <- c("AA", "AH", "AB", "HA", "HH", "HB", "BA", "BH", "BB")
      for(l in gsub("NANA", "", paste0(mR1,mR2))){
        if(l == "AA") observed["AA"] <- observed["AA"] + 1
        if(l == "AH") observed["AH"] <- observed["AH"] + 1
        if(l == "HA") observed["HA"] <- observed["HA"] + 1
        if(l == "AB") observed["AB"] <- observed["AB"] + 1
        
        if(l == "HH") observed["HH"] <- observed["HH"] + 1

        if(l == "BA") observed["BA"] <- observed["BA"] + 1
        if(l == "BH") observed["BH"] <- observed["BH"] + 1
        if(l == "HB") observed["HB"] <- observed["HB"] + 1
        if(l == "BB") observed["BB"] <- observed["BB"] + 1
      }
      lt1 <- table(unlist(mR1))
      lt2 <- table(unlist(mR2))
      expected <- rep(0, 9)
      names(expected) <- c("AA", "AH", "AB", "HA", "HH", "HB", "BA", "BH", "BB")
      for(e1 in c("A", "H", "B")){
        for(e2 in c("A", "H", "B")){
          expected[paste0(e1, e2)] <- (lt1[e1] / sum(lt1)) * (lt2[e2] / sum(lt2))
        }
      }
      expected[is.na(expected)] <- 0
      
      expe <- sum(observed) * expected[sG]
      
      chiSQ <- 0
      for(ijkik in  c("AA", "AB", "BA", "BB")){
        if(expe[ijkik] != 0) chiSQ <- chiSQ + (((observed[ijkik] - expe[ijkik]) ^ 2) / expe[ijkik])
      }
      cat(observed[sG], "-", expe, "chiSQ:", chiSQ, "\n")
      results[y, x] <- chiSQ
      
      
      m1o <-  rep("", length(F2))
      m1o[genotypes[as.character(allRegions[x,"Flanking.Marker"]), F2] == "A"] <- "I/I"
      m1o[genotypes[as.character(allRegions[x,"Flanking.Marker"]), F2] == "H"] <- ""
      m1o[genotypes[as.character(allRegions[x,"Flanking.Marker"]), F2] == "B"] <- "D/D"

      m2o <-  rep("", length(F2))
      m2o[genotypes[as.character(allRegions[y,"Flanking.Marker"]), F2] == "A"] <- "I/I"
      m2o[genotypes[as.character(allRegions[y,"Flanking.Marker"]), F2] == "H"] <- ""
      m2o[genotypes[as.character(allRegions[y,"Flanking.Marker"]), F2] == "B"] <- "D/D"

      m1   <- genotype(m1o)
      m2   <- genotype(m2o)
      tryCatch(
       results2[y, x] <- unlist(LD(m1,m2))$"X^2"
      , error = function(e){ results2[y, x] <- NA; }
      )

    }
  }
}
rownames(results) <- allRegions[,"Flanking.Marker"]
colnames(results) <- allRegions[,"Flanking.Marker"]

bfmipref <- which(allRegions[match(allRegions[,1], colnames(results)),"Prefered.Allele"] == "BFMI")
b6npref <- which(allRegions[match(allRegions[,1], colnames(results)),"Prefered.Allele"] == "B6N")
nopref <- which(allRegions[match(allRegions[,1], colnames(results)),"Prefered.Allele"] == "?")

locations <- seq(0, 1, length.out=ncol(results))

colz <- rep("black", nrow(allRegions))
colz[bfmipref] <- "orange"
colz[b6npref] <- "gray"

namez <- paste0(allRegions[,"Chr"],":", 
          round(allRegions[,"Start"] / 1000000, 0),"-",
          round(allRegions[,"Stop"] / 1000000, 0), " ", 
          allRegions[,"origin"])

nTests <- (ncol(results) * ncol(results)) / 2
op <- par(mar= c(6, 6, 4, 2) + 0.1)
image(-log10(pchisq(results, 1, lower.tail=FALSE)), breaks=c(0, -log10(0.1 / nTests), -log10(0.05 / nTests), -log10(0.01 / nTests), -log10(0.001 / nTests), 100), col=c("white", "gray", "yellow", "orange", "red"), xaxt='n', yaxt='n', main = "Genetic incompatibility (BFMIxB6n)")
axis(1, at=locations[bfmipref], namez[bfmipref], las=2, cex.axis=0.7, col.axis="orange")
axis(1, at=locations[b6npref], namez[b6npref], las=2, cex.axis=0.7, col.axis="gray")
axis(1, at=locations[nopref], namez[nopref], las=2, cex.axis=0.7, col.axis="black")
axis(2, at=locations[bfmipref], namez[bfmipref], las=2, cex.axis=0.7, col.axis="orange")
axis(2, at=locations[b6npref], namez[b6npref], las=2, cex.axis=0.7, col.axis="gray")
axis(2, at=locations[nopref], namez[nopref], las=2, cex.axis=0.7, col.axis="black")
grid(ncol(results),ncol(results))
sx <- 0
for(x in table(allRegions[,"Chr"])){
  l <- (locations[sx] + locations[sx + 1]) / 2
  abline(h=l); abline(v=l)
  sx = sx + x
}
legend("topleft", c("p > 0.1", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01", "p < 0.001"), fill  = c("white", "gray", "yellow", "orange", "red"),cex=0.7, bg="white")
box()

# retrieve entrez gene ID
require("biomaRt")
mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
res <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = c("ensembl_gene_id"), values = unique(as.character(allgenes[,1])), mart = mart)
write.table(res, "genesToEntrez.txt",sep="\t")


# Get all genes in each regions
allgenes <- NULL
allgenes <- read.table("Analysis/Mat_B6NgenesInfo.txt", sep="\t", header=TRUE)
allgenes <- rbind(allgenes, read.table("Analysis/Mat_BFMIgenesInfo.txt", sep="\t", header=TRUE))
allgenes <- rbind(allgenes, read.table("Analysis/Pat_B6NgenesInfo.txt", sep="\t", header=TRUE))
allgenes <- rbind(allgenes, read.table("Analysis/Pat_BFMIgenesInfo.txt", sep="\t", header=TRUE))

# Look at P:P interactions in BIOGRID
dd <- read.csv("stringdb/BIOGRID-ALL-3.4.132.tab2.txt",sep="\t")
dd <- dd[which(dd[,"Organism.Interactor.A"] == 10090 & dd[,"Organism.Interactor.B"] == 10090),]
dd <- dd[which(as.character(dd[,"Entrez.Gene.Interactor.A"]) %in% as.character(res[,2]) & dd[,"Entrez.Gene.Interactor.B"] %in% res[,2]),]
rr <- dd[which(as.character(dd[,c(2,3,8,9)][,3]) != as.character(dd[,c(2,3,8,9)][,4])), c(2,3,8,9)]

rr <- cbind(rr, RegionGeneA = NA, RegionGeneB = NA, ChiSq = NA)

# Which region the gene is in
whichRegion <- function(ensgid){
  qqqq <- allgenes[which(allgenes[,1] == ensgid),]
  gChr <- qqqq[1,"chromosome_name"]
  gLoc <- qqqq[1,"start_position"]
  rownames(allRegions[which(allRegions[,"Chr"] == gChr & as.numeric(allRegions[,"Start"]) < gLoc & as.numeric(allRegions[,"Stop"]) > gLoc ),])
}

## get regions and Chi2 scores
for(x in 1:nrow(rr)){
  ensgidA <- res[which(res[,2] == rr[x,1]),1]
  ensgidB <- res[which(res[,2] == rr[x,2]),1]
  rr[x, "RegionGeneA"] <- whichRegion(ensgidA)[1]
  rr[x, "RegionGeneB"] <- whichRegion(ensgidB)[1]
  if(!is.na(rr[x, "RegionGeneA"]) && !is.na(rr[x, "RegionGeneB"])){
    rr[x, "ChiSq"] <- results[as.numeric(rr[x, "RegionGeneA"]) , as.numeric(rr[x, "RegionGeneB"]) ]
  }
}

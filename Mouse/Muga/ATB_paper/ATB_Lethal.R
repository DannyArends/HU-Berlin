#
# Analyse ATB in megaMuga data from multiple generations after beagle phasing
#

# Load the genotype call data

source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/incompatible.R")

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")

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


regionnames <- paste0(allRegions[,"Chr"],":", 
          round(allRegions[,"Start"] / 1000000, 0),"-",
          round(allRegions[,"Stop"] / 1000000, 0), " ", 
          allRegions[,"origin"])

rownames(allRegions) <- regionnames

## new idea compare 2 regions against each other:
results <- matrix(0, nrow(allRegions), nrow(allRegions))
results2 <- matrix(0, nrow(allRegions), nrow(allRegions))

for(y in 1:nrow(allRegions)){
  mR1 <- genotypes[as.character(allRegions[y,"Flanking.Marker"]), F2]
  chrR1 <- allRegions[y,"Chr"]
  for(x in 1:nrow(allRegions)) {
    mR2 <- genotypes[as.character(allRegions[x,"Flanking.Marker"]), F2]
    chrR2 <- allRegions[x,"Chr"]
    if(chrR1 == chrR2){
      results[x,y] <- NA
      results2[x,y] <- NA
    }else{
      results[x,y] <- incompatible(mR1, mR2, TRUE)$chisq
      #results2[x,y] <- doLD(mR1, mR2)
    }
  }
}
rownames(results) <- regionnames
colnames(results) <- regionnames

write.table(results, "ChiSqScores_Incompatible.txt", sep="\t")

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
results <- read.table("ChiSqScores_Incompatible.txt", sep="\t", check.names=FALSE)
rr <- apply(results,1,as.numeric)
rownames(rr) <- colnames(rr)
results <- rr


# Load the Allele Transmission Biased regions
regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)

allRegions <- rbind(cbind(regionsMat, origin = "M"), cbind(regionsPat, origin = "P"))
allRegions <- allRegions[with(allRegions, order(Chr, Start)), ]

regionnames <- paste0(allRegions[,"Chr"],":", 
          round(allRegions[,"Start"] / 1000000, 0),"-",
          round(allRegions[,"Stop"] / 1000000, 0), " ", 
          allRegions[,"origin"])

rownames(allRegions) <- regionnames

# Generate an overview plot
nTests <- (ncol(results) * ncol(results)) / 2
LODscores <- -log10(pchisq(results, 1, lower.tail=FALSE))
thresholdL <- -log10(0.05 / nTests)
thresholdH <- -log10(0.01 / nTests)

significant <-  names(which(apply(LODscores,1,max,na.rm=TRUE) > thresholdL))
length(significant)
significantH <-  names(which(apply(LODscores,1,max,na.rm=TRUE) > thresholdH))
length(significantH)

allRegions <- allRegions[significant, ]
regionnames <- rownames(allRegions)
results <- results[significant,significant]

allRegions["8:62-63 P", "Prefered.Allele"] <- "B6N"
allRegions["9:88-91 M", "Prefered.Allele"] <- "B6N"
allRegions["11:13-13 M", "Prefered.Allele"] <- "BFMI"

bfmipref <- which(allRegions[,"Prefered.Allele"] == "BFMI")
b6npref <- which(allRegions[,"Prefered.Allele"] == "B6N")
nopref <- which(allRegions[,"Prefered.Allele"] == "?")

colz <- rep("black", nrow(allRegions))
colz[bfmipref] <- "orange"
colz[b6npref] <- "gray"

locations <- seq(0, 1, length.out=ncol(results))

sum(LODscores > -log10(0.1 / nTests), na.rm=TRUE) / 2
sum(LODscores > -log10(0.05 / nTests), na.rm=TRUE) / 2
sum(LODscores > -log10(0.01 / nTests), na.rm=TRUE) / 2

#tiff("genetic_incompatibilities.tiff", width=1024 * 3, height=1024* 3, res=300)
#postscript("Figure_3.eps", horizontal = FALSE, paper = "special", width=16, height=14, pointsize = 18)
  op <- par(mar= c(7, 7, 4, 10) + 0.1, xpd=FALSE)
  image(-log10(pchisq(results, 1, lower.tail=FALSE)), 
        breaks=c(0, -log10(0.1 / nTests), -log10(0.05 / nTests), -log10(0.01 / nTests), -log10(0.001 / nTests), 100), 
        col=c("white",  "darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4"), xaxt='n', yaxt='n', main = "Genetic incompatibilities")

  axis(1, at=locations[bfmipref], regionnames[bfmipref], las=2, cex.axis=1, col.axis="darkorange2")
  axis(1, at=locations[b6npref],  regionnames[b6npref], las=2, cex.axis=1, col.axis="cornflowerblue")
  axis(1, at=locations[nopref],   regionnames[nopref], las=2, cex.axis=1, col.axis="gray")

  axis(2, at=locations[bfmipref], regionnames[bfmipref], las=2, cex.axis=1, col.axis="darkorange2")
  axis(2, at=locations[b6npref],  regionnames[b6npref], las=2, cex.axis=1, col.axis="cornflowerblue")
  axis(2, at=locations[nopref],   regionnames[nopref], las=2, cex.axis=1, col.axis="gray")

  grid(ncol(results),ncol(results))
  sx <- 0
  for(x in table(allRegions[,"Chr"])){
    l <- (locations[sx] + locations[sx + 1]) / 2
    abline(h=l); abline(v=l)
    sx = sx + x
  }
  op <- par(xpd=TRUE)
  legend("topright", c("p > 0.1", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01", "p < 0.001"), 
         fill  = c("white", "darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4"), 
         cex=1, inset=c(-0.18,0))
  legend("topright", c("B6N", "BFMI"), fill = c("cornflowerblue", "orange"), inset=c(-0.18,0.14))
  box()

#dev.off()

# Get all genes in each regions
allgenes <- NULL
allgenes <- read.table("Analysis/Mat_B6NgenesInfo.txt", sep="\t", header=TRUE)
allgenes <- rbind(allgenes, read.table("Analysis/Mat_BFMIgenesInfo.txt", sep="\t", header=TRUE))
allgenes <- rbind(allgenes, read.table("Analysis/Pat_B6NgenesInfo.txt", sep="\t", header=TRUE))
allgenes <- rbind(allgenes, read.table("Analysis/Pat_BFMIgenesInfo.txt", sep="\t", header=TRUE))

# retrieve entrez gene ID
require("biomaRt")
mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
ensembltoentrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = c("ensembl_gene_id"), values = unique(as.character(allgenes[,1])), mart = mart)
write.table(ensembltoentrez, "genesToEntrez.txt",sep="\t")

# Look at P:P interactions in BIOGRID
biogriddata <- read.csv("stringdb/BIOGRID-ALL-3.4.132.tab2.txt",sep="\t")
biogriddata <- biogriddata[which(biogriddata[,"Organism.Interactor.A"] == 10090 & biogriddata[,"Organism.Interactor.B"] == 10090),]
biogriddata <- biogriddata[which(as.character(biogriddata[,"Entrez.Gene.Interactor.A"]) %in% as.character(ensembltoentrez[,2]) & biogriddata[,"Entrez.Gene.Interactor.B"] %in% ensembltoentrez[,2]),]

# Create results
interacting <- biogriddata[which(as.character(biogriddata[,c(2,3,8,9)][,3]) != as.character(biogriddata[,c(2,3,8,9)][,4])), c(2,3,8,9)]
interacting <- cbind(interacting, RegionGeneA = NA, RegionGeneB = NA, ChiSq = NA, LOD = NA)

# Which region the gene is in
geneLocation <- function(allgenes, ensgid, verbose = TRUE){
  ii <- allgenes[which(allgenes[,1] == ensgid),]
  cat(ensgid, ": Found", nrow(ii), "matching genes\n")
  gChr <- ii[1,"chromosome_name"];gLoc <- ii[1,"start_position"]; gLoc1 <- ii[1,"end_position"]  # Chr:Location
  cat(ensgid, ": ", gChr, "-", gLoc, "\n")
  return(paste0(gChr,":", round(gLoc / 1000000,d = 2),"-", round(gLoc1/ 1000000, d = 2)))
}


# Which region the gene is in
whichRegion <- function(allgenes, ensgid, verbose = TRUE){
  ii <- allgenes[which(allgenes[,1] == ensgid),]
  cat(ensgid, ": Found", nrow(ii), "matching genes\n")
  gChr <- ii[1,"chromosome_name"];gLoc <- ii[1,"start_position"]  # Chr:Location
  cat(ensgid, ": ", gChr, "-", gLoc, "\n")
  rownames(allRegions[which(allRegions[,"Chr"] == gChr & as.numeric(allRegions[,"Start"]) < gLoc & as.numeric(allRegions[,"Stop"]) > gLoc ),])
}

## Get regions and Chi2 scores
for(x in 1:nrow(interacting)){
  ensgidA <- ensembltoentrez[which(ensembltoentrez[,2] == interacting[x,1]),1]
  ensgidB <- ensembltoentrez[which(ensembltoentrez[,2] == interacting[x,2]),1]
  interacting[x, "RegionGeneA"] <- whichRegion(allgenes, ensgidA)[1]
  interacting[x, "RegionGeneB"] <- whichRegion(allgenes, ensgidB)[1]
  if(!is.na(interacting[x, "RegionGeneA"]) && !is.na(interacting[x, "RegionGeneB"])){
    interacting[x, "ChiSq"] <- results[ interacting[x, "RegionGeneA"] , interacting[x, "RegionGeneB"] ]
    interacting[x, "LOD"] <- LODscores[ interacting[x, "RegionGeneA"] , interacting[x, "RegionGeneB"] ]
  }
}
colnames(interacting)[1:4] <- c("Entrez.A","Entrez.B","Symbol.A", "Symbol.B")


# Look at P:P interactions in STRINGDB
stringdbdata <- read.csv("stringdb/10090.protein.actions.v10.txt",sep="\t")
ensembltopeptide <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "mgi_symbol"), filters = c("ensembl_gene_id"), values = unique(as.character(allgenes[,1])), mart = mart)
write.table(ensembltopeptide, "genesToPeptide.txt",sep="\t")

stringdbdata <- stringdbdata[as.numeric(stringdbdata[,"score"]) > 700,]
dim(stringdbdata)
stringdbdata <- stringdbdata[which(as.character(stringdbdata[,"item_id_a"]) %in% paste0("10090.",as.character(ensembltopeptide[,2])) & 
                                                stringdbdata[,"item_id_b"] %in% paste0("10090.",as.character(ensembltopeptide[,2]))),]
dim(stringdbdata)
interacting <- stringdbdata[which(stringdbdata[,1] != stringdbdata[,2]),]
interacting <- cbind(interacting, MGI_A = NA, LOC_A = NA, MGI_B = NA, LOC_B = NA, RegionGeneA = NA, RegionGeneB = NA, ChiSq = NA, LOD = NA)

## Get regions and Chi2 scores
for(x in 1:nrow(interacting)){
  ensgidA <- ensembltopeptide[which(paste0("10090.",as.character(ensembltopeptide[,2])) == as.character(interacting[x,1])),1]
  mgiIDA  <- ensembltopeptide[which(paste0("10090.",as.character(ensembltopeptide[,2])) == as.character(interacting[x,1])),3]
  ensgidB <- ensembltopeptide[which(paste0("10090.",as.character(ensembltopeptide[,2])) == as.character(interacting[x,2])),1]
  mgiIDB  <- ensembltopeptide[which(paste0("10090.",as.character(ensembltopeptide[,2])) == as.character(interacting[x,2])),3]
  interacting[x, "RegionGeneA"] <- whichRegion(allgenes, ensgidA)[1]
  interacting[x, "RegionGeneB"] <- whichRegion(allgenes, ensgidB)[1]
  if(!is.na(interacting[x, "RegionGeneA"]) && !is.na(interacting[x, "RegionGeneB"])){
    interacting[x, "ChiSq"] <- results[ interacting[x, "RegionGeneA"] , interacting[x, "RegionGeneB"] ]
    interacting[x, "LOD"] <- LODscores[ interacting[x, "RegionGeneA"] , interacting[x, "RegionGeneB"] ]
    interacting[x, "MGI_A"] <- mgiIDA
    interacting[x, "MGI_B"] <- mgiIDB
    interacting[x, "LOC_A"] <- geneLocation(allgenes, ensgidA)
    interacting[x, "LOC_B"] <- geneLocation(allgenes, ensgidB)
  }
}
write.table(interacting[which(interacting[,"LOD"] > threshold),], "interactionsInRegions_Sign3x3.txt",sep="\t", row.names=FALSE)

significant <- interacting[which(interacting[,"LOD"] > threshold),]


sift[which(paste0("10090.", as.character(sift[, 7])) %in%  c(as.character(interacting[,"item_id_a"]), as.character(interacting[,"item_id_b"]))),]


setwd("E:/Mouse/DNA/Sequencing/BFMI")

if(!file.exists("20140515_VEP_BFMI860mm10.txt")){
  library(openxlsx)
  snpdata <- NULL
  for(x in 1:22){
    snpdata <- rbind(snpdata, read.xlsx("20140515_VEP_BFMI860mm10.xlsx", sheet = x))
  }
  write.table(snpdata, "20140515_VEP_BFMI860mm10.txt", sep = "\t", row.names=FALSE)
}else{
  snpdata <- read.table("20140515_VEP_BFMI860mm10.txt", sep = "\t", header=TRUE)
}

setwd("E:/Mouse/DNA/MegaMuga/")

snpsubset <- snpdata[which(as.character(snpdata[, "Gene_Name"]) %in%  c(as.character(significant[,"MGI_A"]), as.character(significant[,"MGI_B"]))),]
shortlist <- snpsubset[which(grepl("CODING_REGION", snpsubset[,"Region"]) & snpsubset[,"FuncClass"] == "NON_SYNONYMOUS_CODING"), c(1:5,8,12,13,14,15,16,19)]

write.table(shortlist, "mutationsInRegions_Sign3x3.txt",sep="\t", row.names=FALSE)


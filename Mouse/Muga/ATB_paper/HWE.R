#
# Analyse HWE in megaMuga data from multiple generations after beagle phasing
#

source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

# Load the phenotype data
setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                              # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

# Load the genotype call data
setwd("E:/Mouse/DNA/MegaMuga/")
phased.vcf <- read.table(gzfile(paste0("Analysis/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("Analysis/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header

samples <- colnames(phased.vcf)[-c(1:9)]

phased.vcf[, samples] <- apply(phased.vcf[, samples],2,function(x){unlist(lapply(strsplit(x,":"),"[",1))})                   # Only the genotypes
rownames(phased.vcf) <- phased.vcf[,"ID"]
phased.vcf <- phased.vcf[,-which(colnames(phased.vcf) %in% c("ID", "QUAL","FILTER","INFO","FORMAT"))]                              # Remove unneeded columns
phased.vcf[1:10,1:10]

# Change the genotype coding
phased.AHBp <- phased.vcf                                                                         # Copy
phased.AHBp[, samples] <- apply(phased.AHBp[, samples], 2, fromVCFcode.AHBp)                      # Change coding to A H0 H1 B

# Change the genotype coding
phased.geno <- phased.vcf                                                                         # Copy
phased.geno[, samples] <- apply(phased.geno[, samples], 2, fromVCFcode.geno)                      # Change coding to AA AB BA BB

phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% samples),]       # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]

#Load marker annotation
marker.annot <- read.table("Analysis/markerAnnotation.txt", sep="\t")
# Remove the duplicated markers
duplicates <- rownames(marker.annot)[which(duplicated(apply(marker.annot[,c("Chr","Pos")],1,paste0,collapse="")))]
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% duplicates),]

# Order chromosomes in the normal way
chromosomes  <- c(1:19, "X", "Y", "M")
phased.AHBpN <- NULL
phased.genoN <- NULL
for(chr in chromosomes){
  phased.AHBpN <- rbind(phased.AHBpN, phased.AHBp[which(phased.AHBp[,"CHROM"] == chr),])
  phased.genoN <- rbind(phased.genoN, phased.geno[which(phased.geno[,"CHROM"] == chr),])
}
phased.AHBp <- phased.AHBpN
phased.geno <- phased.genoN

generation <- 28
counts <- read.table(paste0("Analysis/TransmissionBias_annotated_",generation,".txt"))

usedmarkers <- rownames(counts)

library(heterozygous)
HWEf2 <- HWE(phased.geno[usedmarkers, F2])

map <- marker.annot[rownames(phased.geno[usedmarkers,F2]), ]
map.autosomes <- map[which(map[,"Chr"] %in% 1:19),]

## Create the chromosome plot
op <- par(mfrow=c(1,1))
ymax <- max(as.numeric(map.autosomes[,"Pos"]))

plot(c(1,19), c(0, ymax), t = 'n', xlab="Chromosome", ylab="Position (Mbp)", xaxt='n', yaxt='n', main="HW equilibrium in generation 28")
axis(1, at=1:19, paste0("Chr ", 1:19), las=2, cex.axis=0.7)
axis(2, at=seq(0,ymax, 20000000), seq(0, ymax, 20000000) / 1000000, las=2, cex.axis=0.9)
for(x in 1:19){
  onChr <- rownames(map.autosomes[which(map.autosomes[,"Chr"] == as.character(x)),])
  colorz <- rep(2, length(HWEf2[onChr]))
  colorz[which(is.nan(HWEf2[onChr]))] <- 1 
  colorz[which(p.adjust(HWEf2[onChr],"BH") < 0.01)] <- 3  
  points(rep(x,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch='-', col=c('gray', 'green', 'red')[colorz], cex=1)
}

HWEdata <- cbind(HWEf2, p.adjust(HWEf2,"BH"))
colnames(HWEdata) <- c("HWP","HWPbh")

write.table(HWEdata, "Analysis/HWEdata28.txt", sep = "\t", quote=FALSE)



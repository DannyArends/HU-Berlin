#
# zcat /home/danny/NAS/Mouse/RNAseqKourosh/AKR-1739/HF3W7DRXY/L001/AKR-1739_S1_L001_R1_001.fastq.gz | head -n 100
# zcat /home/danny/NAS/Mouse/RNAseqKourosh/alignments/trimmed/AKR-1739trimmed.fastq.gz | head -n 100
#

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("BiocParallel")
setwd("/home/danny/NAS/Mouse/RNAseqKourosh/")

fastq <- c(
  "AKR-1739/HF3W7DRXY/L001/AKR-1739_S1_L001_R1_001.fastq.gz",
  "AKR-1749/HF3W7DRXY/L001/AKR-1749_S2_L001_R1_001.fastq.gz",
  "AKR-1756/HF3W7DRXY/L001/AKR-1756_S3_L001_R1_001.fastq.gz",
  "AKR-1761/HF3W7DRXY/L001/AKR-1761_S4_L001_R1_001.fastq.gz",
  
  "B6-1_1009437/HF3W7DRXY/L001/B6-1_1009437_S5_L001_R1_001.fastq.gz",
  "B6-2_1009438/HF3W7DRXY/L001/B6-2_1009438_S6_L001_R1_001.fastq.gz",
  "B6-3_1009442/HF3W7DRXY/L001/B6-3_1009442_S7_L001_R1_001.fastq.gz",
  "B6-4_1009443/HF3W7DRXY/L001/B6-4_1009443_S8_L001_R1_001.fastq.gz",
  
  "BFMI-1_860-29698/HF3W7DRXY/L001/BFMI-1_860-29698_S9_L001_R1_001.fastq.gz",
  "BFMI-2_860-29699/HF3W7DRXY/L001/BFMI-2_860-29699_S10_L001_R1_001.fastq.gz",
  "BFMI-3_860-29700/HF3W7DRXY/L001/BFMI-3_860-29700_S11_L001_R1_001.fastq.gz",
  "BFMI-4_860-29692/HF3W7DRXY/L001/BFMI-4_860-29692_S12_L001_R1_001.fastq.gz",
  
  "SJL-1216/HF3W7DRXY/L001/SJL-1216_S13_L001_R1_001.fastq.gz",
  "SJL-1223/HF3W7DRXY/L001/SJL-1223_S14_L001_R1_001.fastq.gz",
  "SJL-1227/HF3W7DRXY/L001/SJL-1227_S15_L001_R1_001.fastq.gz",
  "SJL-1234/HF3W7DRXY/L001/SJL-1234_S16_L001_R1_001.fastq.gz"
)
infiles <- c()
for(sampleID in 1:length(fastq)){
  fileloc <- "/home/danny/NAS/Mouse/RNAseqKourosh/"
  sampleName <- unlist(lapply(strsplit(fastq[sampleID], "/"), "[",1))
  R1 <- paste0(fileloc, fastq[sampleID])

  outputfolder <- "/home/danny/NAS/Mouse/RNAseqKourosh/alignments/"
  logfile <- paste0(outputfolder, "log", sampleName,".txt")
  referencefolder <- "/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68"
  reference <- paste0(referencefolder, "/GRCm38_68.fa")

  outputBASE <- paste0(outputfolder, sampleName)
  outputSIRBAM <- paste0(outputBASE, ".trimmed.aligned.sorted.dedup.recalibrated.bam")
  infiles <- c(infiles, outputSIRBAM)
}

chrominfo     <- read.table("MouseChrInfo.txt", sep="\t", header=TRUE, colClasses=c("character","integer","logical"))

# Use the whole GFF and then the 4 miRNAs we're intersted in
mouse         <- makeTxDbFromGFF("miRNA.gtf", format = "gtf", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/")
transcriptsByGene <- transcriptsBy(mouse, by = "gene")
bamfiles          <- BamFileList(infiles, yieldSize = 100000)

register(MulticoreParam(workers=8))

se                <- summarizeOverlaps(transcriptsByGene, bamfiles, mode="Union", singleEnd=TRUE, ignore.strand=TRUE)
rawreads <- assay(se)                                                                       # Extract the raw-reads per gene
write.table(rawreads, file="miRNAreads.txt",sep="\t")



### Local
library(preprocessCore)


setwd("D:/Ddrive/Github/HU-Berlin/Kourosh")
allR <- read.table("ALLreads.txt")
miRna <- read.table("miRNAreads.txt")
allD <- rbind(miRna,allR)
nI <- which(apply(allD,1,sum) == 0)
allD <- allD[-nI,]
ix <- allD < 1
allD[ix] <- NA
allDT <- log2(allD + 1)
allDTN <- normalize.quantiles(as.matrix(allDT))
colnames(allDTN) <- gsub(".trimmed.aligned.sorted.dedup.recalibrated.bam", "", colnames(allD))
rownames(allDTN) <- rownames(allD)
boxplot(allDT)
allDTN[ix] <- 0
boxplot(allDTN)

op <- par(mfrow=c(3,1))
for(x in 1:3){
  AKR <- allDTN[x, 1:4]
  B6 <- allDTN[x, 5:8]
  BFMI <- allDTN[x, 9:12]
  SJL <- allDTN[x, 13:16]
  plot(c(0,5), c(0,5), t = 'n', xaxt='n', main = rownames(allDTN)[x], xlab = "Strain", ylab="miRNA expression")
  rect(0.75, 0, 1.25, mean(B6), col = "green")
  points(c(1,1), c(mean(B6) - sd(B6), mean(B6) + sd(B6)), t = 'l')
  points(c(0.85,1.15), c(mean(B6) - sd(B6), mean(B6) - sd(B6)), t = 'l')
  points(c(0.85,1.15), c(mean(B6) + sd(B6), mean(B6) + sd(B6)), t = 'l')

  rect(1.75, 0, 2.25, mean(BFMI), col = "red")
  points(c(2,2), c(mean(BFMI) - sd(BFMI), mean(BFMI) + sd(BFMI)), t = 'l')
  points(c(1.85,2.15), c(mean(BFMI) - sd(BFMI), mean(BFMI) - sd(BFMI)), t = 'l')
  points(c(1.85,2.15), c(mean(BFMI) + sd(BFMI), mean(BFMI) + sd(BFMI)), t = 'l')

  rect(2.75, 0, 3.25, mean(AKR), col = "blue")
  points(c(3,3), c(mean(AKR) - sd(AKR), mean(AKR) + sd(AKR)), t = 'l')
  points(c(2.85,3.15), c(mean(AKR) - sd(AKR), mean(AKR) - sd(AKR)), t = 'l')
  points(c(2.85,3.15), c(mean(AKR) + sd(AKR), mean(AKR) + sd(AKR)), t = 'l')

  rect(3.75, 0, 4.25, mean(SJL), col = "yellow")
  points(c(4,4), c(mean(SJL) - sd(SJL), mean(SJL) + sd(SJL)), t = 'l')
  points(c(3.85,4.15), c(mean(SJL) - sd(SJL), mean(SJL) - sd(SJL)), t = 'l')
  points(c(3.85,4.15), c(mean(SJL) + sd(SJL), mean(SJL) + sd(SJL)), t = 'l')

  axis(1, at = 1:4, c("B6", "BFMI","AKR", "SJL"))
}

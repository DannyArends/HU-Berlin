# RNA Seq - Expression data analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

# Install needed packages form bioconductor.org
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicAlignments", "GenomicFeatures", "Rsamtools"))

# Download GTF file from: ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/
library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")

#setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
chrominfo <- read.table("GTF/MouseChrInfo.txt", sep="\t", header=TRUE, colClasses=c("character","integer","logical"))
sampleIDs <- read.table("FASTQ/sampledescription.txt",sep="\t", header=TRUE)

mouse         <- makeTranscriptDbFromGFF("GTF/Mus_musculus.GRCm38.76.gtf", format = "gtf", exonRankAttributeName="exon_number", 
                                         species="Mus musculus", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/")
exonsByGene   <- exonsBy(mouse, by = "gene")
infiles       <- list.files(pattern="recalibrated.bam$", full=TRUE)
bamfiles      <- BamFileList(infiles, yieldSize = 1000000, asMates=TRUE)
se            <- summarizeOverlaps(exonsByGene, bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

head(assay(se))

namez <- rownames(assay(se))

#RPKM = (10^9 * C)/(N * L)
#C = Number of reads mapped to a gene
#N = Total mapped reads in the experiment
#L = gene length in base-pairs for a gene

exonicGeneSizes <- lapply(exonsByGene, function(x){sum(width(reduce(x)))})
N <- apply(assay(se), 2, sum)

n <- 1
RPKM <- t(apply(assay(se), 1, function(x){
  RPKM <- (10^9 * x) / (N * as.numeric(exonicGeneSizes[n]))
  n <<- n + 1
  return(round(RPKM, d = 3))
}))

corenames <- as.character(sampleIDs[match(as.numeric(unlist(lapply(strsplit(colnames(RPKM),"_"),"[",1))), sampleIDs$Lib_id), "Core.name"])
corenames <- lapply(strsplit(corenames,"_"),"[",1)
colnames(RPKM) <- corenames
rownames(RPKM) <- namez
cat("We called expressions for", nrow(RPKM), "genes\n")

write.table(RPKM, file="RPKM.txt", sep="\t")

RPKM_MPI <- read.csv("MPI_RPKM_ANALYSIS/2014-07-04_BFMI_RPKM.txt", sep="\t", header=TRUE, row.names=1)       # RNA-Seq summary data from MDC
cat("MPI called expressions for", nrow(RPKM_MPI), "genes\n")
RPKM_MPI <- RPKM_MPI[, unlist(lapply(colnames(RPKM), grep, colnames(RPKM_MPI)))]

RPKM_MPI <- RPKM_MPI[which(rownames(RPKM_MPI) %in% rownames(RPKM)),]
RPKM_GAB <- RPKM[which(rownames(RPKM) %in% rownames(RPKM_MPI)),]

png("comparison_RPKM_GAB_MPI.png")
  op <- par(mfrow=c(1,2))
  plot(RPKM_GAB[,1], RPKM_MPI[,1], main= colnames(RPKM_GAB)[1], log = "xy")
  plot(RPKM_GAB[,2], RPKM_MPI[,2], main= colnames(RPKM_GAB)[2], log = "xy")
dev.off()

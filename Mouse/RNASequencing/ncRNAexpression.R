# RNA Seq - Expression data analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
# Chromosome info structure, holding the number of chromosomes, their names and the length per chromosome
chrominfo   <- read.table("GTF/MouseChrInfo.txt", sep="\t", header=TRUE, colClasses=c("character","integer","logical"))
geneinfo    <- read.table("GTF/Mus_musculus.GRCm38.76.gtf", sep="\t", header=TRUE)
allmiRNA    <- unique(unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(geneinfo[geneinfo[,2] == "miRNA",9]),";"),"[",1))," "),"[",2)))

RPKM <- read.table("Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")
RPKM <- RPKM[which(RPKM[,"Mean.B6NxBFMI860.12.L"] > 1 | RPKM[,"Mean.BFMI860.12xB6N.L"] > 1),]
RPKM <- RPKM[which(RPKM[,"tTest_F1"] < 0.1),]

write.table(RPKM[which(RPKM[,"ensembl_gene_id"] %in% allmiRNA),], "Analysis/miRNA_RPKM_Qnorm_ANN_AddDom.txt", sep = "\t")

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
# Chromosome info structure, holding the number of chromosomes, their names and the length per chromosome
chrominfo     <- read.table("GTF/MouseChrInfo.txt", sep="\t", header=TRUE, colClasses=c("character","integer","logical"))

# Download GTF file from: ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/
mouse         <- makeTranscriptDbFromGFF("GTF/Mus_musculus.GRCm38.76.gtf", format = "gtf", exonRankAttributeName="exon_number", 
                                         species="Mus musculus", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/")

transcriptsBy <- exonsBy(mouse, by = "exon")
infiles       <- list.files(path = "./Analysis", pattern="recalibrated.bam$", full=TRUE)
bamfiles      <- BamFileList(infiles, yieldSize = 1000000, asMates=TRUE)
se            <- summarizeOverlaps(exonsByGene, bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

head(assay(se))                                                                             # Show the first 10 lines of the data matrix
rawreads <- assay(se)                                                                       # Extract the raw-reads per gene

write.table(rawreads, file="Analysis/RawReadsPerExon.txt", sep="\t")                               # Save the raw reads to a file

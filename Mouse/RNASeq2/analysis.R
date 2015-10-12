# Pipeline for paired end RNA-Seq analysis
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Okt, 2015
# first written Okt, 2015

# Install needed packages form bioconductor.org (only needs to be done once)
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GenomicAlignments", "GenomicFeatures", "Rsamtools","preprocessCore"))

baminputdir     <- as.character("/home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/FASTQ/")
reference.gtf   <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.81.gtf"                      # wget ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz
short.gtf       <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.81.short.gtf"
sampledescrpath <- "/home/arends/RNASeq/Experiment/SampleDescription.txt"
chrdescrpath    <- "/home/arends/Mouse/GTF/MouseChrInfo.txt"

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")

chrs <- c(1:19,"X","Y","MT")
sampledescr  <- read.table(sampledescrpath, sep="\t", header=TRUE, colClasses="character")
chrominfo    <- read.table(chrdescrpath, sep="\t", header=TRUE, colClasses=c("character","integer","logical"))
chrominfo    <- chrominfo[chrominfo[,1] %in% chrs,]

baminputfolder <- dir(baminputdir)
inputbfiles <- baminputfolder[which(grepl("output", baminputfolder))]
sampleNames <- substr(inputbfiles, 1, 4)
cat("Found ", length(sampleNames), " samples in folder: ", baminputdir, "\n")

bamfiles <- NULL
for(msample in inputbfiles){
  files <- dir(paste0(baminputdir, msample))
  fname <- files[which(grepl(".aln.sort.rmdup.rg.realigned.recal.bam", files))]
  if(!((length(fname) == 0) && (typeof(fname) == "character"))) bamfiles <- c(bamfiles, paste0(baminputdir, msample, "/", fname))
}                                                                                                             #done <- substr(bamfiles, 1, 4)
cat("Found ", length(bamfiles), " aligned bam files in folder: ", baminputdir, "\n")

if(!file.exists(short.gtf)){
  gtf <- read.table(reference.gtf,sep="\t")
  cat(readLines(reference.gtf,n=5),sep="\n", file=short.gtf)
  write.table(gtf[gtf[,1] %in% chrs,], file=short.gtf,sep="\t", row.names = FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
}

if(!file.exists("/home/share/genomes/mm10/exonsByGene.RData")){
  mouseTDB    <- makeTranscriptDbFromGFF(short.gtf, format = "gtf", exonRankAttributeName="exon_number", 
                                         species="Mus musculus", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/")
  save(mouseTDB, file = "/home/share/genomes/mm10/mouseTDB.RData")
}else{
  load("/home/share/genomes/mm10/mouseTDB.RData")
}

if(!file.exists("/home/share/genomes/mm10/exonsByGene.RData")){
  exonsByGene <- exonsBy(mouseTDB, by = "gene")
  save(exonsByGene, file = "/home/share/genomes/mm10/exonsByGene.RData")
}else{
  load("/home/share/genomes/mm10/exonsByGene.RData")
}

library(BiocParallel)
register(MulticoreParam(workers = 2))

bamfilelist <- BamFileList(bamfiles, yieldSize = 1000000, asMates = TRUE)
sT <- proc.time()
overlap     <- summarizeOverlaps(features = exonsByGene, reads = bamfilelist, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
(proc.time() - sT)[3]

write.table(assay(overlap), file="Reads18Samples.txt", sep="\t")


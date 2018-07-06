
library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")

setwd("/home/danny/RNAseqBFMItranscripts/")

# Chromosome info structure, holding the number of chromosomes, their names and the length per chromosome
reference.gtf <- "/home/danny/genomes/mm10/Mus_musculus.GRCm38.89.gtf"
chrominfo     <- read.table("/home/danny/genomes/mm10/MouseChrInfo.txt", sep="\t", header=TRUE, colClasses=c("character","integer","logical"))

# Download GTF file from: ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/
mouse         <- makeTxDbFromGFF(reference.gtf, format = "gtf", 
                                         organism="Mus musculus", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-89/gtf/mus_musculus/")

exonsByTranscripts <- exonsBy(mouse, by = "tx", use.names = TRUE)

dirs <- list.files(path = "/home/danny/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/FASTQ/", pattern=".output", full=TRUE)
files <- NULL
for(d in dirs){
  files <- c(files, list.files(path = d, pattern="recal.bam$", full=TRUE))
}

bamfiles      <- BamFileList(files, yieldSize = 1000000, asMates=TRUE)
se            <- summarizeOverlaps(exonsByTranscripts, bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

head(assay(se))                                                                             # Show the first 10 lines of the data matrix
rawreads <- assay(se)                                                                       # Extract the raw-reads per gene

write.table(rawreads, file="transcript_expression.txt", sep="\t")                           # Save the raw reads to a file


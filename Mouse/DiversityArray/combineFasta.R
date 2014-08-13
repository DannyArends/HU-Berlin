# Analysis of JAX mouse diversity chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written June, 2014

setwd("E:/Mouse/DNA/")

chromosomes <- c(1:19, "X", "Y", "MT")

# Create the database fasta from the ENSEMBLE fasta files
cat("", file="DiversityArray/Analysis/Mus_musculus.GRCm38.74.dna.fasta")
for(chr in chromosomes){
  fastadata <- readLines(paste0("Annotation/GenomeSequence/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa"))
  cat(fastadata, sep="\n", file="DiversityArray/Analysis/Mus_musculus.GRCm38.74.dna.fasta", append = TRUE)
}

# After this create the BLAST database files:
# makeblastdb -in Mus_musculus.GRCm38.74.dna.fasta -dbtype nucl -title Mus_musculus.GRCm38 -out Mus_musculus.GRCm38.74.dna.db

setwd("E:/Mouse/DNA/DiversityArray/")

# Create the blast query fasta from the JAX mouse diversity chip annotation files
chrAnnotationJAX <- NULL
aa <- lapply(chromosomes, function(chr){
  chrAnnotation <- read.table(paste0("Annotation/ProbeAnnotation/chr",chr,".txt"), header=FALSE, sep='\t')        # SNP / Chromosome annotation of the mouse diversity CHIP
  chrAnnotationJAX <<- rbind(chrAnnotationJAX, chrAnnotation[,c(1, 9)])                                           # Take only the annotation of interest
})
colnames(chrAnnotationJAX) <- c("JAX_ID", "Sequence")

emptySequence <- which(nchar(as.character(chrAnnotationJAX[,"Sequence"])) == 0)                     # Remove the one without a sequence
chrAnnotationJAX <- chrAnnotationJAX[-emptySequence, ]

dupEntries <- which(duplicated(as.character(chrAnnotationJAX[,"JAX_ID"])))                          # Remove the duplicate entries
chrAnnotationJAX <- chrAnnotationJAX[-dupEntries,]

cat("", file="Analysis/JAXsequences.fasta")
apply(chrAnnotationJAX, 1, function(annotation){
  cat(paste0(">", annotation["JAX_ID"], "\n", as.character(annotation["Sequence"]), "\n"), file="Analysis/JAXsequences.fasta", append = TRUE)
})

# After this use blastn to query the database for the location of the JAX probe sequences
# blastn -task blastn -query JAXsequences.fasta -db Mus_musculus.GRCm38.74.dna.db -perc_identity 100 -outfmt 6 -evalue=0.1 -out JAXblasted.txt

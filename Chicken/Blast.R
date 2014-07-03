# Analysis of Chicken 600K SNP chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

setwd("E:/Chicken/DNA/")
chromosomes <- c(1:28, 32, "W", "Z", "MT")

# Create the database fasta from the ENSEMBLE fasta files
cat("", file="600KSNPChip/Analysis/Gallus_gallus.Galgal4.74.dna.fasta")
for(chr in chromosomes){
  fastadata <- readLines(paste0("Annotation/GenomeSequence/Gallus_gallus.Galgal4.74.dna.chromosome.", chr, ".fa"))
  cat(fastadata, sep="\n", file="600KSNPChip/Analysis/Gallus_gallus.Galgal4.74.dna.fasta", append = TRUE)
  cat("Done with chromosome", chr, "\n")
}

# After this create the BLAST database files:
# makeblastdb -in Gallus_gallus.Galgal4.74.dna.fasta -dbtype nucl -title Gallus_gallus.Galgal4 -out Gallus_gallus.Galgal4.74.dna.db

setwd("E:/Chicken/DNA/600KSNPChip/")

arrayAnnotation <- read.table("Annotation/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", header=TRUE, colClasses="character", na.string="---")
snpAlleles <- as.character(apply(apply(arrayAnnotation[, c("Allele.A","Allele.B")], 1, sort), 2, paste0, collapse="/"))
snpAlleles <- paste("[", snpAlleles, "]",sep="")

nProbes <- nrow(arrayAnnotation)
cat("", file="Analysis/ArraySequences.fasta")
for(x in 1:nProbes){
  dnasequence <- gsub(snpAlleles[x],"N", arrayAnnotation[x, "Flank"],fixed = TRUE)
  cat(paste0(">", arrayAnnotation[x,"Probe.Set.ID"], "\n", dnasequence, "\n"), file="Analysis/ArraySequences.fasta",append=TRUE)
  if(x %% 1000 == 0)cat("Done", x,"/",nProbes, "\n")
}

# After this use blastn to query the database for the location of the JAX probe sequences
# blastn -task blastn -query ArraySequences.fasta -db Gallus_gallus.Galgal4.74.dna.db -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 5 -out ProbeLocationBLAST.txt

# Create a full fasta file for the genome build GalGal4 (ensembl version 78)
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015
#

chromosomes <- c(1:28,32, "Z", "W", "MT")

setwd("E:/Chicken/DNA")

cat("", file = "Reference/Gallus_gallus.Galgal4.fasta")
for(chr in chromosomes){
  cat(readLines(paste0("Reference/Gallus_gallus.Galgal4.dna.chromosome.",chr,".fa")), sep="\n", file = "Reference/Gallus_gallus.Galgal4.fasta", append=TRUE)
}

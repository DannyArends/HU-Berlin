setwd("D:/Edrive/Goat/DNA/Sequencing")
library(seqinr)
library(Biostrings)

fastaSeq <- read.fasta("genome/GoatGenome_6_8_14.fa")

setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
allsnps <- read.table(file="SNPs.filtered.annotated.txt", sep = "\t", header=TRUE)
snps <- read.table(file="SNPsInCDS.filtered.annotated.txt", sep = "\t", header=TRUE)

toIUPAC <- function(ref, alt){
  sorted <- sort(c(ref,alt))
  if(sorted[1] == "A"){
    if(sorted[2] == "G") return("R")
    if(sorted[2] == "C") return("M")
    if(sorted[2] == "T") return("W")
  }
  if(sorted[1] == "C"){
    if(sorted[2] == "G") return("S")
    if(sorted[2] == "T") return("Y")
  }
  if(sorted[1] == "G"){
    if(sorted[2] == "T") return("K")
  }
}

for(x in 1:nrow(snps)){
  chr <- snps[x,"CHROM"]
  pos <- snps[x,"POS"]
  bpsbefore <- (pos-60):(pos-1)
  bpsafter <- (pos+1):(pos+60)
  seqbefore <- fastaSeq[[chr]][bpsbefore]
  snpsBefore <- which(allsnps[,"CHROM"] == chr & allsnps[,"POS"] %in% bpsbefore)
  if(length(snpsBefore) > 0){
    for(y in snpsBefore){
      pos <- allsnps[y,"POS"]
      iupac <- toIUPAC(as.character(allsnps[y,"REF"]), as.character(allsnps[y,"ALT"]))
      seqbefore[which(bpsbefore == pos)] <- iupac
    }
  }
  seqafter <- fastaSeq[[chr]][bpsafter]
  snpsAfter <- which(allsnps[,"CHROM"] == chr & allsnps[,"POS"] %in% bpsafter)
  if(length(snpsAfter) > 0){
    for(y in snpsAfter){
      pos <- allsnps[y,"POS"]
      iupac <- toIUPAC(as.character(allsnps[y,"REF"]), as.character(allsnps[y,"ALT"]))
      seqafter[which(bpsafter == pos)] <- iupac
    }
  }
  snps[x, "Sequence"] <- toupper(paste0(paste0(seqbefore, collapse=""), "[", snps[x, "REF"],"/", snps[x, "ALT"], "]", paste0(seqafter, collapse="")))
}

write.table(snps, "SNPsSequenceKASP_IUPAC.txt", sep='\t', quote = FALSE, row.names=FALSE)

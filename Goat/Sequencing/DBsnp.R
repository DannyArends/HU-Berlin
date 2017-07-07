#
# DBSNP annotation of SNPs found using capture Seq in Capra Hircus
# DB snp uses the CHIR 1.0 genome
#


setwd("D:/Edrive/Goat/DNA/Sequencing")
library(seqinr)
library(Biostrings)

fastaSeq <- read.fasta("genome/GoatGenome_6_8_14.fa")

## get the fas files from DBSNP
## Concatenate them together: type rs_ch6.fas rs_ch8.fas rs_ch14.fas > dbsnp_3chr.fasta
## Make the blasy database
## makeblastdb -in dbsnp_3chr.fasta -dbtype nucl -title dbsnp_3chr -out dbsnp_3chr.db

setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
snps <- read.csv("filteredSNPs.txt", comment.char="#", sep="\t", header=TRUE, check.names=FALSE)

toIUPAC <- function(ref, alt){
iupac <- rbind( c("R", "A", "G"),
                c("Y", "C", "T"),
                c("S", "G", "C"),
                c("W", "A", "T"),
                c("K", "G", "T"),
                c("M", "A", "C"))
  A <- ref == iupac[,2] | ref == iupac[,3]
  B <- alt == iupac[,2] | alt == iupac[,3]
  return(iupac[A & B, 1])
}

get51BP <- function(chr, pos, ref, alt){
  front <- fastaSeq[[chr]][(pos - 26) :(pos-1)]
  base <- toIUPAC(ref,alt)
  back <- fastaSeq[[chr]][(pos + 1) : (pos+26)]
  return(paste0(paste0(front,collapse=""),base, paste0(back,collapse="")))
}

cat("", file="inputSeq.fasta")
for(x in 1:nrow(snps)){
  cat(">", paste0(as.character(snps[x,"CHROM"]),".", snps[x,"POS"]), "\n", file="inputSeq.fasta", append=TRUE)
  cat(get51BP(snps[x,"CHROM"],snps[x,"POS"],snps[x,"REF"],snps[x,"ALT"]), "\n", file="inputSeq.fasta", append=TRUE)
}

#Blast the input against the DBSNP file:
#blastn -task blastn -query inputSeq.fasta -db ../DBsnp/dbsnp_3chr.db -perc_identity 99 -outfmt 6 -evalue 0.1 -num_alignments 5 -out locations.txt

blastResults <- read.table("locations.txt", colClasses="character")
blastResults[-which(blastResults[,"V4"] != 51),]                        # Matches with less then 51 basepairs

for(x in 1:nrow(snps)){
  snpID <- paste0(as.character(snps[x,"CHROM"]),".", snps[x,"POS"])
  blastR <- which(blastResults[,1] == snpID)
  if(length(blastR) >= 1){
     dbSNPids <- gsub("gnl|dbSNP|","",blastResults[blastR,2], fixed=TRUE)
     snps[x,"ID"] <- paste0(dbSNPids, collapse=";")
  }
}

write.table(snps, "filteredSNPs_DBsnp.txt", sep="\t")
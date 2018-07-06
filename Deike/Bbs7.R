library(seqinr)
library(Biostrings)

setwd("D:/Edrive/Mouse/DNA/")

fastaSeq <- read.fasta("Annotation/Mus_musculus.GRCm38.dna.chromosome.3.fa.gz")
gff3file <- gzfile("Annotation/Mus_musculus.GRCm38.92.chromosome.3.gff3.gz","r")
gff3 <- read.csv(gff3file, comment.char="#", sep="\t", header=FALSE)
close(gff3file)
gff3exons <- gff3[which(gff3[, "V3"] == "exon"), ]                  # Get only the rows containing exons from the GFF3 file
snps <- read.csv("Annotation/Bbs7_region_output.raw.snps.indels.all.vcf", comment.char="#", sep="\t", header=TRUE, check.names=FALSE, colClasses="character")


exons <- c("ENSMUSE00001286926", "ENSMUSE00001222852", "ENSMUSE00000172582", "ENSMUSE00000395912", "ENSMUSE00000368075")
exons <- gff3exons[unlist(lapply(lapply(exons, grep,  gff3exons[,9]), "[", 1)), ]
exons <- exons[sort(exons[,4],index.return=TRUE)$ix,]

sequences <- vector("list", 5)

for(x in 1:nrow(exons)){
  exonsequence <- fastaSeq[[1]][exons[x, 4]:exons[x, 5]]
  names(exonsequence) <- as.character(exons[x, 4]:exons[x, 5])
  sequences[[x]] <- toupper(exonsequence)
}

IUPAC = rbind(R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"), W = c("A", "T"), K = c("G", "T"), M = c("A", "C"))

for(x in 1:length(sequences)){
  SeqisSNP <- names(sequences[[x]])[which(names(sequences[[x]]) %in% as.character(snps[,"POS"]))]
  if(length(SeqisSNP) > 0){
    for(pos in SeqisSNP){
      cat("exons", 6-x, " snp at", pos, "\n")
      inSNPs <- which(snps[, "POS"] == pos)
      snpcode <- sort(as.character(c(snps[inSNPs, "REF"], snps[inSNPs, "ALT"])))
      iupaccode <- rownames(IUPAC)[which(IUPAC[,1] == snpcode[1] & IUPAC[,2] == snpcode[2])]
      cat("snpcode", snpcode[1], ":", snpcode[2], "=", iupaccode, "\n")
      sequences[[x]][SeqisSNP] <- iupaccode
    }
  }
}
sequences[[3]][1:4] <- rep("X",4)

write.table(rbind(names(unlist(sequences)), unlist(sequences)), "SequenceEx5toEx1_bbs7.txt", col.names=FALSE,sep="\t", row.names=FALSE)


setwd("D:/Edrive/Goat/DNA/Sequencing")
library(seqinr)
library(Biostrings)

fastaSeq <- read.fasta("genome/GoatGenome_6_8_14.fa")
gff3 <- read.csv("Genome/ugff3.gff3", comment.char="#", sep="\t", header=FALSE)

setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
snps <- read.csv("filteredSNPs.txt", comment.char="#", sep="\t", header=TRUE, check.names=FALSE)

genes <- c("CSN1S1", "CSN2", "CSN1S2", "LOC102178810", "CSN3", "STC1", "HSF1", "DGAT1")

genedata <- vector("list", length(genes))
names(genedata) <- genes

for(gene in genes){
  gffall <- gff3[grep(gene, gff3[,"V9"]),]
  gffexons <- gffall[which(gffall[,"V3"] == "exon"),]
  gffcds <- gffall[which(gffall[,"V3"] == "CDS"),]
  proteins <- unique(unlist(lapply(strsplit(as.character(gffcds[,"V9"]),";"), function(x){
    gsub("protein_id=", "", x[grep("protein_id=",x)])
  })))
  
  genedata[[gene]] <- vector("list", length(proteins))
  names(genedata[[gene]]) <- proteins
  for(protein in proteins) {
    gffproteincds <- gffcds[grep(protein, gffcds[,"V9"]),]
    cat(gene, protein, nrow(gffproteincds), "\n")
    gffproteincds <- cbind(gffproteincds, sequence = NA)
    gffproteincds <- gffproteincds[sort(gffproteincds[,"V4"], index.return=TRUE)$ix,]
    for(x in 1:nrow(gffproteincds)){
      gffproteincds[x, "sequence"] = paste0(fastaSeq[[gffproteincds[x, "V1"]]][gffproteincds[x, "V4"]:gffproteincds[x, "V5"]], collapse="")
    }
    genedata[[gene]][[protein]] <- gffproteincds
  }
}

toProtein <- function(proteingff){
  codingsequence <- c()
  for(x in 1:nrow(proteingff)){
    codingsequence <- paste0(codingsequence, proteingff[x, "sequence"], collapse="")
  }
  if(all(proteingff[,"V7"] == "+")){
    return(list(codingsequence, as.character(translate(DNAString(codingsequence)))))
  }else{
    return(list(codingsequence, as.character(translate(reverseComplement(DNAString(codingsequence))))))
  }
}

toProtein(gendata[[1]][[1]])
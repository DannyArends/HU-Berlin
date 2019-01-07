#
# Annotation of SNPs found using capture Seq in Capra Hircus
# NOTE: This approach doesnt work since it assumes all CDS regions are alligned on a 3 codon DNA barrier
#
setwd("D:/Edrive/Goat/DNA/Sequencing")
library(seqinr)
library(Biostrings)

fastaSeq <- read.fasta("genome/GoatGenome_6_8_14.fa")
gff3 <- read.csv("Genome/ugff3.gff3", comment.char="#", sep="\t", header=FALSE)
gff3genes <- gff3[which(gff3[, "V3"] == "gene"), ]                  # Get only the rows containing genes from the GFF3 file
gff3exons <- gff3[which(gff3[, "V3"] == "exon"), ]                  # Get only the rows containing exons from the GFF3 file
gff3cds <- gff3[which(gff3[, "V3"] == "CDS"), ]                     # Get only the rows containing CDS from the GFF3 file

setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
snps <- read.csv("filteredSNPs.txt", comment.char="#", sep="\t", header=TRUE, check.names=FALSE)

snps <- cbind(snps, inGene = NA)                                    # Add a column to hold the answer for each SNP
snps <- cbind(snps, geneOri = NA)                                   # Add a column to hold the answer for each SNP
snps <- cbind(snps, locGene = NA)                                   # Add a column to hold the answer for each SNP
for (x in 1:nrow(snps)) {                                           # Go through each of the SNPs
  chr <- snps[x,"CHROM"]                                            # Chromosome the SNP is located on
  pos <- as.numeric(snps[x,"POS"])                                  # Position of the SNP

  # Test if the SNP is within the start and end position of any of the GFF genes
  inR <- which(gff3genes[, "V1"] == chr & pos > as.numeric(gff3genes[,"V4"]) & pos < as.numeric(gff3genes[,"V5"]))

  # inR will contain the row number of the gene in the GFF file (no gene: it will be NULL)
  if(length(inR) > 0) {    
    msplit <- strsplit(as.character(gff3genes[inR, 9]), ";")[[1]]
    geneName <- gsub("Name=", "", msplit[grep("Name=", msplit)])
    snps[x, "inGene"] <- geneName                                             # We did find a gene around our SNP
    snps[x, "geneOri"] <- as.character(gff3genes[inR, 7])                     # Orientation
    if(snps[x, "geneOri"] == "+"){
      snps[x, "locGene"] <- 1 + (pos - as.numeric(gff3genes[inR, "V4"]))        # Location (pos strand)
    }else{
      snps[x, "locGene"] <- 1 + (as.numeric(gff3genes[inR, "V5"]) - pos)        # Location (neg strand)
    }
  } else {
    snps[x, "inGene"] <- NA # We did not find a gene around our SNP
  }
}

snps <- cbind(snps, inExon = NA)                # Add a column to hold the answer for each snp
snps <- cbind(snps, locExon = NA)               # Add a column to hold the answer for each snp
for (x in 1:nrow(snps)) {                       # Go through each of the SNPs
  if(!is.na(snps[x, "inGene"])) {               # if we are not in a gene, we cannot be in an exon
    chr <- snps[x,"CHROM"]                      # Chromosome the SNP is located on
    pos <- as.numeric(snps[x,"POS"])            # Position of the SNP

    # Test if the SNP is within the start and end position of any of the GFF exons
    inR <- which(gff3exons[, "V1"] == chr & pos > as.numeric(gff3exons[,"V4"]) & pos < as.numeric(gff3exons[,"V5"]))
    # inR will contain the row number of the exon in the GFF file (no gene: it will be NULL)
    if(length(inR) > 0) {    
      snps[x, "inExon"] <- 1   # We did find an exon around our SNP
      if(snps[x, "geneOri"] == "+"){
        snps[x, "locExon"] <- 1 + (pos - as.numeric(gff3exons[inR[1], "V4"]))
      }else{
        snps[x, "locExon"] <- 1 + (as.numeric(gff3exons[inR[1], "V5"]) - pos)
      }
    } else {
      snps[x, "inExon"] <- NA  # We did not find an exon around our SNP
    }
  }else{
    snps[x, "inExon"] <- NA  # Not in a gene not in an exon
  }
}

toComplement <- function(x){
  complement <- c()
  for(e in x){
    if(e == "t") complement <- c(complement, "a")
    if(e == "a") complement <- c(complement, "t")
    if(e == "c") complement <- c(complement, "g")
    if(e == "g") complement <- c(complement, "c")
  }
  return(complement)
}

AAtranslate <- rbind(c("G", "Glycine", "Gly"), c("P", "Proline", "Pro"),
      c("A", "Alanine", "Ala"), c("V", "Valine", "Val"),
      c("L", "Leucine", "Leu"), c("I", "Isoleucine", "Ile"),
      c("M", "Methionine", "Met"), c("C", "Cysteine", "Cys"),
      c("F", "Phenylalanine", "Phe"),c("Y", "Tyrosine", "Tyr"),
      c("W", "Tryptophan", "Trp"), c("H", "Histidine", "His"),
      c("K", "Lysine", "Lys"), c("R", "Arginine", "Arg"),
      c("Q", "Glutamine", "Gln"), c("N", "Asparagine", "Asn"),
      c("E", "Glutamic Acid", "Glu"),c("D", "Aspartic Acid", "Asp"),
      c("S", "Serine", "Ser"),c("T", "Threonine", "Thr"))

rownames(AAtranslate) <- AAtranslate[,1]
AAtranslate <- AAtranslate[,-1]

snps <- cbind(snps, inCDS = NA)          # Add a column to hold the answer for each snp
snps <- cbind(snps, locCDS = NA)         # Add a column to hold the answer for each snp
snps <- cbind(snps, locAA = NA)          # Add a column to hold the answer for each snp
snps <- cbind(snps, seqAA = NA)          # Add a column to hold the answer for each snp
snps <- cbind(snps, refAA = NA)          # Add a column to hold the answer for each snp
snps <- cbind(snps, newAA = NA)          # Add a column to hold the answer for each snp
snps <- cbind(snps, altAA = NA)          # Add a column to hold the answer for each snp
for (x in 1:nrow(snps)) {                   # Go through each of the SNPs
  if(!is.na(snps[x, "inGene"])) {             # if we are not in a gene, we cannot be in an exon
    chr <- snps[x,"CHROM"]                      # Chromosome the SNP is located on
    pos <- as.numeric(snps[x,"POS"])            # Position of the SNP

    # Test if the SNP is within the start and end position of any of the GFF exons
    inR <- which(gff3cds[, "V1"] == chr & pos > as.numeric(gff3cds[,"V4"]) & pos < as.numeric(gff3cds[,"V5"]))
    # inR will contain the row number of the exon in the GFF file (no gene: it will be NULL)
    if(length(inR) > 0) {
      snps[x, "inCDS"] <- 1   # We did find an exon around our SNP
      if(snps[x, "geneOri"] == "+") {
        snps[x, "locCDS"] <- 1 + (pos - as.numeric(gff3cds[inR[1], "V4"]))
      } else {
        snps[x, "locCDS"] <- 1 + (as.numeric(gff3cds[inR[1], "V5"]) - pos)
      }
      snps[x, "locAA"] <- snps[x, "locCDS"] %% 3
      if(snps[x, "locAA"] == 0) snps[x, "locAA"] = 3
      
      if(snps[x, "geneOri"] == "+") {
        codonS <- pos - (snps[x, "locAA"] - 1)
        codon <- fastaSeq[[chr]][codonS:(codonS + 2)]
        snps[x, "seqAA"] <- toupper(paste0(codon, collapse=""))
        codon[snps[x, "locAA"]] <- tolower(snps[x, "ALT"])
        snps[x, "newAA"] <- toupper(paste0(codon,collapse=""))
      }else{
        codonS <- pos + (snps[x, "locAA"] - 1)
        codon <-fastaSeq[[chr]][codonS:(codonS - 2)]
        snps[x, "seqAA"] <- toupper(paste0(toComplement(codon),collapse=""))
        codon[snps[x, "locAA"]] <- tolower(snps[x, "ALT"])
        snps[x, "newAA"] <- toupper(paste0(toComplement(codon),collapse=""))
      }
      refAA <- as.character(translate(DNAString(snps[x, "seqAA"])))
      newAA <- as.character(translate(DNAString(snps[x, "newAA"])))
      snps[x, "refAA"] <- paste0(AAtranslate[refAA, 1], " (", refAA, ")")
      snps[x, "altAA"] <- paste0(AAtranslate[newAA, 1], " (", newAA, ")")
    } else {
      snps[x, "inCDS"] <- NA  # We did not find an exon around our SNP
    }
  }else{
    snps[x, "inCDS"] <- NA  # Not in a gene not in an exon
  }
}

snps[,c(1, 2, 4, 5, 43:54)]

snpsincds <- snps[which(snps[,"inCDS"] == 1),c(1, 2, 4, 5, 43:54)]
write.table(snpsincds, file="SNPsInCDS.filtered.annotated.txt", sep = "\t", quote = FALSE, row.names=FALSE)
write.table(snps[,c(1, 2, 4, 5, 43:54)], file="SNPs.filtered.annotated.txt", sep = "\t", quote = FALSE, row.names=FALSE)

snps <- cbind(snps, MAF = NA)          # Add a column to hold the answer for each snp
for(x in 1:nrow(snps)){
  tbl <- table(unlist(strsplit(as.character(unlist(snps[x, 10:42])), "")))
  snps[x, "MAF"] = round(((tbl / sum(tbl)) * 100)[as.character(snps[x, "ALT"])], 1)
}

# K23 == K33 
samples <- read.table("../samples_description.txt",sep="\t", row.names = 1)

rownames(samples)[!rownames(samples) %in% colnames(snps)]
colnames(snps[,10:42])[!colnames(snps[,10:42]) %in% rownames(samples)]

for(breed in unique(samples[,1])){
  ind <- rownames(samples)[which(samples[,1] == breed)]
  snps <- cbind(snps, NA)
  colnames(snps)[ncol(snps)] <- paste0("MAF", breed)
  for(x in 1:nrow(snps)){
    tbl <- table(unlist(strsplit(as.character(unlist(snps[x, ind])), "")))
    snps[x, paste0("MAF", breed)] = ((tbl / sum(tbl)) * 100)[as.character(snps[x, "ALT"])]
    if(is.na(snps[x, paste0("MAF", breed)])) snps[x, paste0("MAF", breed)] = 0
    snps[x, paste0("MAF", breed)] <- round(as.numeric(snps[x, paste0("MAF", breed)]), 1)
  }
}

snps[1:5,]
write.table(snps, file="SNPs.filtered.annotated.mafs.txt", sep = "\t", quote = FALSE, row.names=FALSE)
write.table(snps[,-c(3,7,8,9)], file="SNPs.filtered.annotated.mafs.noinfo.txt", sep = "\t", quote = FALSE, row.names=FALSE,na="")

genes <- unique(na.omit(snps[,"inGene"]))
exonIDs <- NULL
for(gene in genes){
  allexons <- gff3exons[grep(gene, gff3exons[,"V9"]), ]#c(1,2,3,4,5,7)]
  allexons <- allexons[sort(allexons[,"V4"], index.return=TRUE)$ix,]
  allexons <- unique(allexons)
  if(any("BestRefSeq" %in% allexons[,"V2"])){
    allexons <- allexons[allexons[,"V2"] == "BestRefSeq",]
  }else{
    # No refseq, just take the transcript X1
    ii <- grep("X1", allexons[,"V9"])
    if(length(ii) > 1){
      allexons <- allexons[ii,]
    }
  }
  if(all(allexons[,"V7"] == "+")){
    allexons <- cbind(allexons, Gene = gene, exonN = 1:nrow(allexons))
  }else{
    allexons <- cbind(allexons, Gene = gene, exonN = nrow(allexons):1)
  }
  exonIDs <- rbind(exonIDs, allexons)
}

for (x in 1:nrow(snps)) {                   # Go through each of the SNPs
  if(!is.na(snps[x, "inExon"])){
    chr <- snps[x,"CHROM"]                      # Chromosome the SNP is located on
    pos <- as.numeric(snps[x,"POS"])            # Position of the SNP

    inR <- which(exonIDs[, "V1"] == chr & pos > as.numeric(exonIDs[,"V4"]) & pos < as.numeric(exonIDs[,"V5"]))
    if(length(inR) == 0) {
      cat(x, "is not in a reference exon\n")
      snps[x, "inExon"] <- "?"
    }else if(length(inR) == 1) {
      cat(x, "is in 1 exon: ",inR, exonIDs[inR, "exonN"],"\n")
      snps[x, "inExon"] <- exonIDs[inR, "exonN"]
    }else{
      cat(x, "is in 2 exons\n")
    }
  }
}
idxS <- which(colnames(snps) %in% rownames(samples))

write.table(snps[,-c(3, 7, 8, 9)], file="SNPs.filtered.annotated.mafs.noinfo.exonID.txt", sep = "\t", quote = FALSE, row.names=FALSE,na="")
write.table(snps[,-c(3, 7, 8, 9, idxS)], file="SNPs.filtered.annotated.mafs.noinfo.nogeno.exonID.txt", sep = "\t", quote = FALSE, row.names=FALSE,na="")

#
# Annotation of SNPs found using capture Seq in Capra Hircus
# Using the full CDS to determine which AA is hit, this does not assume DNA codons are aligned with the start of the CDS
# GFF3 from : ftp://ftp.ncbi.nlm.nih.gov/genomes/Capra_hircus/GFF/ref_ASM170441v1_top_level.gff3.gz
# GoatGenome_6_8_14 ASSEMBLY NAME: ARS1

setwd("D:/Edrive/Goat/DNA/Sequencing")
library(seqinr)
library(Biostrings)

fastaSeq <- read.fasta("genome/GoatGenome_6_8_14.fa")
gff3 <- read.csv("Genome/ugff3.gff3", comment.char="#", sep="\t", header=FALSE)

setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
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
  for(protein in proteins) {    # Proteins in this respect is equal to transcripts
    gffproteincds <- gffcds[grep(protein, gffcds[,"V9"]),]
    cat(gene, protein, nrow(gffproteincds), "\n")
    gffproteincds <- cbind(gffproteincds, sequence = NA)
    gffproteincds <- gffproteincds[sort(gffproteincds[,"V4"], index.return=TRUE)$ix,]     # Sort them based on increasing start position
    for(x in 1:nrow(gffproteincds)) {                                                     # Use the genomic fasta to get the positive sequence of the CDS
      gffproteincds[x, "sequence"] = paste0(fastaSeq[[gffproteincds[x, "V1"]]][gffproteincds[x, "V4"]:gffproteincds[x, "V5"]], collapse="")
    }
    genedata[[gene]][[protein]] <- gffproteincds
  }
}

# Get the protein sequence
toProtein <- function(proteingff){
  codingsequence <- c()
  for(x in 1:nrow(proteingff)){
    codingsequence <- paste0(codingsequence, proteingff[x, "sequence"], collapse="")
  }
  if(all(proteingff[,"V7"] == "+")) { # Positive strand gene, just translate the CDS
    return(list(codingsequence, as.character(translate(DNAString(codingsequence), if.fuzzy.codon="solve"))))
  }else{                              # Negative strand gene, reverse complement the DNA and translate the rev complemented CDS
    return(list(tolower(as.character(reverseComplement(DNAString(codingsequence)))), as.character(translate(reverseComplement(DNAString(codingsequence)), if.fuzzy.codon="solve"))))
  }
}

# Used to translate AA codes to pretty names
AAtranslate <- rbind(c("G", "Glycine", "Gly"), c("P", "Proline", "Pro"),
      c("A", "Alanine", "Ala"), c("V", "Valine", "Val"),
      c("L", "Leucine", "Leu"), c("I", "Isoleucine", "Ile"),
      c("M", "Methionine", "Met"), c("C", "Cysteine", "Cys"),
      c("F", "Phenylalanine", "Phe"),c("Y", "Tyrosine", "Tyr"),
      c("W", "Tryptophan", "Trp"), c("H", "Histidine", "His"),
      c("K", "Lysine", "Lys"), c("R", "Arginine", "Arg"),
      c("Q", "Glutamine", "Gln"), c("N", "Asparagine", "Asn"),
      c("E", "Glutamic Acid", "Glu"),c("D", "Aspartic Acid", "Asp"),
      c("S", "Serine", "Ser"),c("T", "Threonine", "Thr"),c("*", "Stop", "Stop"))
rownames(AAtranslate) <- AAtranslate[,1]
AAtranslate <- AAtranslate[,-1]


#toProtein(genedata[[1]][[1]])

snps <- read.csv("filteredSNPs_DBsnp.txt", comment.char="#", sep="\t", header=TRUE, check.names=FALSE)
snps <- cbind(snps, inGene = NA)            # Name of the gene in which the SNP is located
snps <- cbind(snps, GeneOri = NA)           # Orientation of the gene
snps <- cbind(snps, inTranscript = NA)      # Transcript of the gene in which the SNP is located
snps <- cbind(snps, inCDS = NA)             # Coding Sequence in which the gene is located
snps <- cbind(snps, locCDS = NA)            # Location of the SNP in the coding sequence in which the SNP is located
snps <- cbind(snps, locCDSf = NA)           # Location of the SNP in the full coding sequence
snps <- cbind(snps, locAA = NA)             # Location fo the AA codon being hit 1, 2, 3
snps <- cbind(snps, AApos = NA)             # Amino acid position in which the SNP is located
snps <- cbind(snps, Synonymous = NA)        # Does the SNP cause an AA code change ?
snps <- cbind(snps, seqAA = NA)             # Reference Codon
snps <- cbind(snps, refAA = NA)             # Reference AA
snps <- cbind(snps, newAA = NA)             # Codon with SNP in it
snps <- cbind(snps, altAA = NA)             # Alternative AA

SNPsinCDS <- NULL
for(x in 1:nrow(snps)){
  chr <- snps[x, "CHROM"]
  pos <- as.numeric(snps[x,"POS"])                                  # Position of the SNP

  for(gene in names(genedata)){                                     # Go through the genes 1 by 1
    for(transcript in names(genedata[[gene]])){                     # Go through the transcripts of a gene 1 by 1
      tdata <- genedata[[gene]][[transcript]]
      for(cds in 1:nrow(tdata)){
        if(tdata[cds,"V1"] == chr && pos > tdata[cds,"V4"] && pos < tdata[cds,"V5"]){ # >= or <=
          cat("SNP", x, "in", gene, transcript, cds, "\n")
          SNPsinCDS <- rbind(SNPsinCDS, snps[x, ])
          posInTable <- nrow(SNPsinCDS)
          
          # Store: Gene, Orientation, Transcript
          SNPsinCDS[posInTable, "inGene"] <- gene
          SNPsinCDS[posInTable, "GeneOri"] <- as.character(unique(tdata[,"V7"]))
          SNPsinCDS[posInTable, "inTranscript"] <- transcript
          
          # If the gene is on the negative strand the number of the CDS is counted inverse
          if(as.character(unique(tdata[,"V7"])) == "+") {
            SNPsinCDS[posInTable, "inCDS"] <- paste0(cds, " / ", nrow(tdata))
          }else{
            SNPsinCDS[posInTable, "inCDS"] <- paste0((nrow(tdata) - cds) + 1, " / ", nrow(tdata))
          }
          
          # If the gene is on the negative strand the position is counted from the length of the CDS
          if(as.character(unique(tdata[,"V7"])) == "+") {
            SNPsinCDS[posInTable, "locCDS"]<- 1 + (pos - as.numeric(tdata[cds,"V4"]))
          }else{
            SNPsinCDS[posInTable, "locCDS"]<- 1 + (as.numeric(tdata[cds,"V5"]) - pos)
          }
          
          # Figure out how long the previous CDS-es were, again this depends on the orientation of the gene
          if(as.character(unique(tdata[,"V7"])) == "+") {
            if(cds > 1){
              prevCDS <- nchar(paste0(tdata[1:(cds-1), "sequence"], collapse=""))
            }else{
              prevCDS <- 0
            }
          }else{
            if(cds < nrow(tdata)){
              prevCDS <- nchar(paste0(tdata[(cds+1):nrow(tdata), "sequence"], collapse=""))
            }else{
              prevCDS <- 0
            }
          }

          # Get the DNA and protein sequence for this transcript
          proteinData <- toProtein(tdata)

          SNPsinCDS[posInTable, "locCDSf"] <- prevCDS + SNPsinCDS[nrow(SNPsinCDS), "locCDS"]

          # position of the codon of the Amino acid being changed
          SNPsinCDS[posInTable, "locAA"] <- SNPsinCDS[nrow(SNPsinCDS), "locCDSf"] %% 3
          if(SNPsinCDS[posInTable, "locAA"] == 0) SNPsinCDS[nrow(SNPsinCDS), "locAA"] = 3

          # Amino acid position int he protein that is changed
          SNPsinCDS[posInTable, "AApos"] <- ceiling(SNPsinCDS[nrow(SNPsinCDS), "locCDSf"] / 3)

          # Stop and start positions for the amino acid codon, toProtein does the reverseComplement of negative genes
          aaStop <- as.numeric(SNPsinCDS[posInTable, "AApos"]) * 3
          aaStart <- aaStop - 2
          
          # Update the AA position to include the total number of AA's in this transcript
          SNPsinCDS[posInTable, "AApos"] <- paste0(SNPsinCDS[posInTable, "AApos"], " / ", length(strsplit(proteinData[[2]][1], "")[[1]]))
          
          # Get the reference codon hit by the SNP
          codon <- strsplit(proteinData[[1]][1], "")[[1]][aaStart:aaStop]
          codon[SNPsinCDS[posInTable, "locAA"]] <- toupper(codon[SNPsinCDS[posInTable, "locAA"]])
          SNPsinCDS[posInTable, "seqAA"] <- paste0(codon, collapse="")
          refAA <- as.character(translate(DNAString(SNPsinCDS[posInTable, "seqAA"])))
          SNPsinCDS[posInTable, "refAA"] <- paste0(AAtranslate[refAA, 1], " (", refAA, ")")

          # Update the codon based on the SNP (or it's complement)
          if(as.character(unique(tdata[,"V7"])) == "+") {
            codon[SNPsinCDS[posInTable, "locAA"]] <- toupper(SNPsinCDS[nrow(SNPsinCDS), "ALT"])
          }else{
            codon[SNPsinCDS[posInTable, "locAA"]] <- as.character(complement(DNAString(tolower(SNPsinCDS[nrow(SNPsinCDS), "ALT"]))))
          }
          SNPsinCDS[posInTable, "newAA"] <- paste0(codon, collapse="")
          altAA <- as.character(translate(DNAString(SNPsinCDS[posInTable, "newAA"])))
          SNPsinCDS[posInTable, "altAA"] <- paste0(AAtranslate[altAA, 1], " (", altAA, ")")

          # Is this mutation (non-)synonymous ?
          if(SNPsinCDS[posInTable, "refAA"] == SNPsinCDS[posInTable, "altAA"]){
            SNPsinCDS[posInTable, "Synonymous"] <- "Yes"
          }else{
            SNPsinCDS[posInTable, "Synonymous"] <- "Non"
          }
        }
      }
    }
  }
}

SNPIDs <- paste(SNPsinCDS[,"CHROM"], SNPsinCDS[, "POS"],sep="-")

# Write out the table of SNPs in Coding regions
write.table(SNPsinCDS[,-c(3,6,7,8,9,10:42)], file="SNPsInTranscripts.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Some basic statistics on SNPs in coding sequences
table(SNPsinCDS[,"locAA"])
table(SNPsinCDS[,"Synonymous"])


## Annotation of MAFs of the different SNPs in populations
samples <- read.table("../samples_description.txt",sep="\t", row.names = 1) # K23 == K33 

rownames(samples)[!rownames(samples) %in% colnames(SNPsinCDS)]
colnames(SNPsinCDS[,10:42])[!colnames(SNPsinCDS[,10:42]) %in% rownames(samples)]

for(breed in unique(samples[,1])){
  ind <- rownames(samples)[which(samples[,1] == breed)]
  SNPsinCDS <- cbind(SNPsinCDS, NA)
  colname <- paste0("MAF", breed)
  colnames(SNPsinCDS)[ncol(SNPsinCDS)] <- colname
  for(x in 1:nrow(SNPsinCDS)){
    tbl <- table(unlist(strsplit(as.character(unlist(SNPsinCDS[x, ind])), "")))
    SNPsinCDS[x, colname] = ((tbl / sum(tbl)) * 100)[as.character(SNPsinCDS[x, "ALT"])]
    if(is.na(SNPsinCDS[x, paste0("MAF", breed)])) SNPsinCDS[x, colname] = 0
    SNPsinCDS[x, colname] <- round(as.numeric(SNPsinCDS[x, colname]), 1)
  }
}

write.table(SNPsinCDS[,-c(6,7,8,9)], file="SNPsInTranscriptsMAFGeno.txt", sep="\t", quote=FALSE, row.names=FALSE, na = "")
write.table(SNPsinCDS[,-c(6,7,8,9,10:42)], file="SNPsInTranscriptsMAF.txt", sep="\t", quote=FALSE, row.names=FALSE, na = "")


SNPselection <- rbind(
  c("SNP1", "CSN1S1", "ENA|CM004567|CM004567.1", 85981710, "[C/A]"),
  c("SNP1a", "CSN1S1", "ENA|CM004567|CM004567.1", 85988705, "[G/A]"),
  c("SNP2", "CSN2", "ENA|CM004567|CM004567.1", 86008016, "[A/G]"),
  c("SNP3", "CSN1S2", "ENA|CM004567|CM004567.1", 86079098, "[T/C]"),
  c("SNP4", "CSN3", "ENA|CM004567|CM004567.1", 86208960, "[A/G]"),
  c("SNP5", "DGAT1", "ENA|CM004575|CM004575.1", 81331681, "[C/A]"),
  c("SNP6", "HSF1", "ENA|CM004575|CM004575.1", 81326608, "[T/G]")
)

SNPselection <- cbind(SNPselection, reference = NA)
SNPselection <- cbind(SNPselection, front = NA)
SNPselection <- cbind(SNPselection, SNP = NA)
SNPselection <- cbind(SNPselection, back = NA)

for(x in 1:nrow(SNPselection)){
  SNPselection[x, "reference"] <- toupper(paste0(fastaSeq[[ SNPselection[x,3] ]][(as.numeric(SNPselection[x,4]) - 60):(as.numeric(SNPselection[x,4]) + 60)], collapse=""))
  SNPselection[x, "front"] <- toupper(paste0(fastaSeq[[ SNPselection[x,3] ]][(as.numeric(SNPselection[x,4]) - 60):(as.numeric(SNPselection[x,4]) -1)], collapse=""))
  SNPselection[x, "back"] <- toupper(paste0(fastaSeq[[ SNPselection[x,3] ]][(as.numeric(SNPselection[x,4]) + 1):(as.numeric(SNPselection[x,4]) + 60)], collapse=""))
  SNPselection[x, "SNP"] <- SNPselection[x,5]
}

write.table(SNPselection, "FlankingSequences.txt", sep ="\t", row.names=FALSE)


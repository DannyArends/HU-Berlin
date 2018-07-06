setwd("~/genomes/EquCab2")
library(seqinr)
library(Biostrings)
library(biomaRt) 

fastaSeq <- read.fasta("Equus_caballus.EquCab2.dna.fa.gz")
gff3 <- read.csv(gzfile("Equus_caballus.EquCab2.90.gff3.gz"), comment.char="#", sep="\t", header=FALSE)
setwd("/home/danny/RNAseqHorses/")
mdata <- read.table("GenesOfInterest.txt", sep='\t', header=TRUE)

ensemblIDS <- as.character(mdata[,"Ensemble"])
bio.mart <- useMart("ensembl", dataset="ecaballus_gene_ensembl")

rownames(mdata) <- ensemblIDS

res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name"), 
                        filters = c("ensembl_gene_id"), values = ensemblIDS, mart = bio.mart)

rownames(res.biomart) <- res.biomart[,1]
mall <- cbind(mdata, res.biomart[rownames(mdata),])

genes <- unique(mall[,"external_gene_name"])
genes <- genes[-which(genes == "")]

genedata <- vector("list", length(genes))
names(genedata) <- genes
for(gene in genes){
  gffall <- gff3[grep(gene, gff3[,"V9"]),]
  gffmrna <- gffall[which(gffall[,"V3"] == "mRNA"),]
  
  proteins <- unique(unlist(lapply(strsplit(as.character(gffmrna[,"V9"]),";"), function(x){
    gsub("transcript_id=", "", x[grep("transcript_id=",x)])
  })))
  
  genedata[[gene]] <- vector("list", length(proteins))
  names(genedata[[gene]]) <- proteins
  for(protein in proteins) {    # Proteins in this respect is equal to transcripts
    gffall <- gff3[grep(protein, gff3[,"V9"]),]
    gffcds <- gffall[which(gffall[,"V3"] == "CDS"),]
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

SNPIDs <- paste(SNPsinCDS[,"CHROM"], SNPsinCDS[, "POS"],sep="-")

colnames(SNPsinCDS)[1:17] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "11", "14", "15", "16", "21", "22", "23", "24")

# Write out the table of SNPs in Coding regions
write.table(SNPsinCDS, file="SNPsInTranscripts_all.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(SNPsinCDS[,-c(3,7,8,9)], file="SNPsInTranscripts_test.txt", sep="\t", quote=FALSE, row.names=FALSE)

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


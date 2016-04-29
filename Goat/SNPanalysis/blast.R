# Preprocessing of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016


# 2 genome builds (CHIR1.0 2011 RefSeq == GenBank & CHIR2.0 2015 GenBank)
#makeblastdb -in GCF_000317765.1_CHIR_1.0_genomic.fna -dbtype nucl -title GCF_000317765.1_CHIR_1.0 -out GCF_000317765.1_CHIR_1.0_genomic.db
#makeblastdb -in GCA_000317765.2_CHIR_2.0_genomic.fna -dbtype nucl -title GCA_000317765.2_CHIR_2.0 -out GCA_000317765.2_CHIR_2.0_genomic.db

setwd("E:/Goat/DNA/SihamRawData")
ILMNannot <- read.csv("goatiggc_cons_60k_11589869_A.csv", skip=7, colClasses="character")
ILMNannot <- cbind(ILMNannot, Sequence = gsub("\\[[A-Z]/[A-Z]\\]", "N", ILMNannot[,"TopGenomicSeq"]))

if(!file.exists("probes.fasta")){
  cat("", file="probes.fasta")
  for(x in 1:nrow(ILMNannot)){
    if(ILMNannot[x,"Sequence"] != ""){
      cat(">", ILMNannot[x,"Name"], "\n", sep="", file="probes.fasta",append=TRUE)
      cat(as.character(ILMNannot[x,"Sequence"]), "\n", sep="", file="probes.fasta",append=TRUE)
    }
  }
}

# Blast against the 2 reference genomes
#blastn -task blastn -query SihamRawData/probes.fasta -db annotation/CHIR1.0/GCF_000317765.1_CHIR_1.0_genomic.db -perc_identity 95 -outfmt "6 qseqid sseqid sstart send sstrand" -evalue 0.1 -num_alignments 5 -out SihamAnalysis/ProbeLocationCHIR1.0.txt
#blastn -task blastn -query SihamRawData/probes.fasta -db annotation/CHIR2.0/GCA_000317765.2_CHIR_2.0_genomic.db -perc_identity 95 -outfmt "6 qseqid sseqid sstart send sstrand" -evalue 0.1 -num_alignments 5 -out SihamAnalysis/ProbeLocationCHIR2.0.txt

### CHIR 1.0 (2011)
setwd("E:/Goat/DNA/")
toCHR <- read.csv("annotation/CHIR1.0/toCHR.txt", sep = "\t",header=FALSE)

chir1 <- read.csv("SihamAnalysis/ProbeLocationCHIR1.0.txt", sep="\t")
colnames(chir1) <- c("decoder_ID", "Chr", "Start", "Stop", "Strand")                                      # Name the columns
chir1 <- cbind(chir1, length = abs(chir1[,"Start"] - chir1[,"Stop"]))
chir1 <- chir1[which(chir1[,"length"] == 120),]                                                           # Only keep the full matches to the reference

multimapping <- as.character(chir1[which(duplicated(chir1[,"decoder_ID"])),"decoder_ID"])                 # Multiple 100 % matches to the reference
if(length(multimapping) > 0) chir1 <- chir1[-which(chir1[,"decoder_ID"] %in% multimapping),]
chir1 <- cbind(chir1, chrN = toCHR[match(chir1[,"Chr"], toCHR[,2]), 1])                                   # Translate chromosome names from NCXXX -> 1,2,3,X,MT
chir1 <- chir1[-which(is.na(chir1[,"chrN"])), ]                                                           # When a probe is not on a chromosome we remove it
write.table(chir1, "SihamAnalysis/FilteredLocationCHIR1.0.txt", sep="\t",row.names=FALSE,quote=FALSE)

### CHIR 2.0 (2015)
setwd("E:/Goat/DNA/")
toCHR <- read.csv("annotation/CHIR2.0/toCHR.txt", sep = "\t",header=FALSE)

chir2 <- read.csv("SihamAnalysis/ProbeLocationCHIR2.0.txt", sep="\t")
colnames(chir2) <- c("decoder_ID", "Chr", "Start", "Stop", "Strand")                                      # Name the columns
chir2 <- cbind(chir2, length = abs(chir2[,"Start"] - chir2[,"Stop"]))
chir2 <- chir2[which(chir2[,"length"] == 120),]                                                           # Only keep the full matches to the reference

multimapping <- as.character(chir2[which(duplicated(chir2[,"decoder_ID"])),"decoder_ID"])                 # Multiple 100 % matches to the reference
if(length(multimapping) > 0) chir2 <- chir2[-which(chir2[,"decoder_ID"] %in% multimapping),]
chir2 <- cbind(chir2, chrN = toCHR[match(chir2[,"Chr"], toCHR[,2]), 1])                                   # Translate chromosome names from NCXXX -> 1,2,3,X,MT
chir2 <- chir2[-which(is.na(chir2[,"chrN"])), ]                                                           # When a probe is not on a chromosome we remove it
write.table(chir2, "SihamAnalysis/FilteredLocationCHIR2.0.txt", sep="\t",row.names=FALSE,quote=FALSE)


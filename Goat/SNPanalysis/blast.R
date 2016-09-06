# Preprocessing of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016

# Genome build LWLT01: http://www.ncbi.nlm.nih.gov/nuccore/LWLT00000000.1
# Assembly: http://www.ncbi.nlm.nih.gov/assembly/GCA_001704415.1/

# 3 genome builds (CHIR1.0 2011 RefSeq == GenBank, CHIR2.0 2015 GenBank, LWLT01 (San Clemente))
#makeblastdb -in GCF_000317765.1_CHIR_1.0_genomic.fna -dbtype nucl -title GCF_000317765.1_CHIR_1.0 -out GCF_000317765.1_CHIR_1.0_genomic.db
#makeblastdb -in GCA_000317765.2_CHIR_2.0_genomic.fna -dbtype nucl -title GCA_000317765.2_CHIR_2.0 -out GCA_000317765.2_CHIR_2.0_genomic.db
#makeblastdb -in GCA_001704415.1_ASM170441v1_genomic.fna -dbtype nucl -title GCA_001704415.1_ASM170441v1 -out GCA_001704415.1_ASM170441v1_genomic.db

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

# Blast against the 3 reference genomes
#blastn -task blastn -query SihamRawData/probes.fasta -db annotation/CHIR1.0/GCF_000317765.1_CHIR_1.0_genomic.db -perc_identity 95 -outfmt "6 qseqid sseqid sstart send sstrand" -evalue 0.1 -num_alignments 5 -out SihamAnalysis/ProbeLocationCHIR1.0.txt
#blastn -task blastn -query SihamRawData/probes.fasta -db annotation/CHIR2.0/GCA_000317765.2_CHIR_2.0_genomic.db -perc_identity 95 -outfmt "6 qseqid sseqid sstart send sstrand" -evalue 0.1 -num_alignments 5 -out SihamAnalysis/ProbeLocationCHIR2.0.txt
#blastn -task blastn -query SihamRawData/probes.fasta -db annotation/LWLT01/GCA_001704415.1_ASM170441v1_genomic.db -perc_identity 95 -outfmt "6 qseqid sseqid sstart send sstrand" -evalue 0.1 -num_alignments 5 -out SihamAnalysis/ProbeLocationLWLT01.txt

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

### LWLT01 (2016)
setwd("E:/Goat/DNA/")
toCHR <- read.csv("annotation/LWLT01/toCHR.txt", sep = "\t",header=FALSE)

lwlt01 <- read.csv("SihamAnalysis/ProbeLocationLWLT01.txt", sep="\t")
colnames(lwlt01) <- c("decoder_ID", "Chr", "Start", "Stop", "Strand")                                      # Name the columns
lwlt01 <- cbind(lwlt01, length = abs(lwlt01[,"Start"] - lwlt01[,"Stop"]))
lwlt01 <- lwlt01[which(lwlt01[,"length"] == 120),]                                                         # Only keep the full matches to the reference

multimapping <- as.character(lwlt01[which(duplicated(lwlt01[,"decoder_ID"])),"decoder_ID"])                # Multiple 100 % matches to the reference
if(length(multimapping) > 0) lwlt01 <- lwlt01[-which(lwlt01[,"decoder_ID"] %in% multimapping),]
lwlt01 <- cbind(lwlt01, chrN = toCHR[match(lwlt01[,"Chr"], toCHR[,2]), 1])                                 # Translate chromosome names from NCXXX -> 1,2,3,X,MT

unplaced <- which(is.na(lwlt01[,"chrN"]))
if(length(unplaced) > 0) lwlt01 <- lwlt01[-unplaced, ]                                                     # When a probe is not on a chromosome we remove it
write.table(lwlt01, "SihamAnalysis/FilteredLocationLWLT01.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Combine:
setwd("E:/Goat/DNA/")
chir1  <- read.table("SihamAnalysis/FilteredLocationCHIR1.0.txt", sep="\t", row.names=1, header=TRUE)
chir2  <- read.table("SihamAnalysis/FilteredLocationCHIR2.0.txt", sep="\t", row.names=1, header=TRUE)
lwlt01 <- read.table("SihamAnalysis/FilteredLocationLWLT01.txt", sep="\t", row.names=1, header=TRUE)

decoders <- c(rownames(chir1), rownames(chir2), rownames(lwlt01))
combined <- names(which(table(decoders) == 3))

locs.orig <- cbind(chir1[combined, c("chrN", "Start", "Stop", "Strand")], 
              chir2[combined, c("chrN", "Start", "Stop", "Strand")], 
              lwlt01[combined, c("chrN", "Start", "Stop", "Strand")])
colnames(locs.orig) <- c("chir1_Chr", "chir1_Start", "chir1_Stop", "chir1_Strand",
                    "chir2_Chr", "chir2_Start", "chir2_Stop", "chir2_Strand",
                    "lwlt01_Chr", "lwlt01_Start", "lwlt01_Stop", "lwlt01_Strand")

cat("Probes suitable for all 3 genomes:", nrow(locs.orig), "\n")

invertLWLT01Chr <- function(locs, chrs){
  for(chr in chrs){
    onChr <- rownames(locs[which(as.character(locs[,"lwlt01_Chr"]) == chr),])
    maxL <- max(max(locs[onChr, "lwlt01_Start"]),max(locs[onChr, "lwlt01_Stop"]))
    newStart <- maxL - locs[onChr, "lwlt01_Start"]
    newStop <- newStart +  (locs[onChr, "lwlt01_Start"] -  locs[onChr, "lwlt01_Stop"])
    locs[onChr, "lwlt01_Start"] <- newStart
    locs[onChr, "lwlt01_Stop"] <- newStop
  }
  return(locs)
}
locs <- invertLWLT01Chr(locs.orig, c("2", "3", "4", "7", "10", "12", "14", "15", "17", "23", "26", "27", "28"))

locs <- locs[with(locs, order(chir1_Chr, chir1_Start)), ]

getOverlap <- function(on_g1, on_g2) {
  mgroup <- on_g1[which(on_g1 %in% on_g2)]
  return(mgroup)
}

getDiff <- function(on_g1, on_g2, pos_g1, pos_g2){
  mgroup <- getOverlap(on_g1, on_g2)
  dif <- pos_g1[mgroup] - pos_g2[mgroup]
  return(dif)
}

hasBigDiff <- function(on_g1, on_g2, pos_g1, pos_g2){
  mgroup <- getOverlap(on_g1, on_g2)
  dif <- getDiff(on_g1, on_g2, pos_g1, pos_g2)
  return(names(which(abs(diff(dif)) > 10000000)))
  
#  return(mgroup[which(dif > (med + 4 * std) | dif < (-med - 4 * std))])
}

signDiffs <- function(on_g1, on_g2, pos_g1, pos_g2){
  return(diff(sign(getDiff(on_g1, on_g2, pos_g1, pos_g2))))
}

hasFlips <- function(on_g1, on_g2, pos_g1, pos_g2){
  mgroup <- getOverlap(on_g1, on_g2)
  signdiffs <- which(diff(sign(getDiff(on_g1, on_g2, pos_g1, pos_g2))) != 0)
  idx <- as.numeric(unlist(lapply(signdiffs,function(x){return(c(x, x + 1))})))
  return(mgroup[idx])
}

chrs.all <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
          "15", "16", "17", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29")

plotID <- 1
#for(plotID in seq(1, 30, 10)){
  #pdf(paste0("plot", plotID,".pdf"), paper="special", height=8, width=14)
  
  chrs <- chrs.all[plotID:min((plotID+9),length(chrs.all))]
  maxbp <- max(max(apply(locs[,grepl("Start", colnames(locs))],2,max)), max(apply(locs[,grepl("Stop", colnames(locs))],2,max)))
  chrXaxt <- 2
  plot(c(0, (length(chrs) * 4) + 1), c(0, maxbp), t = 'n', xaxt='n', yaxt='n', xlab="Chromosome", ylab = "Mbp")
  for(chr in chrs){
    on_c1 <- rownames(locs[which(as.character(locs[,"chir1_Chr"]) == chr),])
    on_c2 <- rownames(locs[which(as.character(locs[,"chir2_Chr"]) == chr),])
    on_l1 <- rownames(locs[which(as.character(locs[,"lwlt01_Chr"]) == chr),])
    
    pos_c1 <- (locs[on_c1,"chir1_Start"] + locs[on_c1,"chir1_Stop"]) / 2
    names(pos_c1) <- on_c1
    pos_c2 <- (locs[on_c2,"chir2_Start"] + locs[on_c2,"chir2_Stop"]) / 2
    names(pos_c2) <- on_c2
    pos_l1 <- (locs[on_l1,"lwlt01_Start"] + locs[on_l1,"lwlt01_Stop"]) / 2
    names(pos_l1) <- on_l1

    points(x = c(chrXaxt-1, chrXaxt-1), y = c(min(pos_c1), max(pos_c1)), lwd=5, t='l')
    points(x = c(chrXaxt+0, chrXaxt-0), y = c(min(pos_c2), max(pos_c2)), lwd=5, t='l')
    points(x = c(chrXaxt+1, chrXaxt+1), y = c(min(pos_l1), max(pos_l1)), lwd=5, t='l')

    bigD12 <- hasFlips(on_c1, on_c2, pos_c1, pos_c2)
    difD12 <- getDiff(on_c1, on_c2, pos_c1, pos_c2)
    bigD23 <- hasFlips(on_c2, on_l1, pos_c2, pos_l1)
    bigD13 <- hasFlips(on_c1, on_l1, pos_c1, pos_l1)

    for(probe in unique(c(on_c1, on_c2, on_l1))) {
      atC1 <- which(on_c1 == probe)
      atC2 <- which(on_c2 == probe)
      atL1 <- which(on_l1 == probe)
      
      if(length(atC1) > 0 && length(atC2) > 0){
        colz <- as.numeric(probe %in% bigD12) + 1
        if(colz != 1) points(x = c(chrXaxt-1, chrXaxt+0), y = c(pos_c1[atC1], pos_c2[atC2]), lwd=colz, t='l', lty = c("dotted", "solid")[colz], col=c("orange", "red")[colz])
      }
      if(length(atC2) > 0 && length(atL1) > 0){
        colz <- as.numeric(probe %in% bigD23) + 1
        if(colz != 1) points(x = c(chrXaxt+0, chrXaxt+1), y = c(pos_c2[atC2], pos_l1[atL1]), lwd=colz, t='l', lty = c("dotted", "solid")[colz], col=c("purple", "red")[colz])
      }
#      if(length(atC1) > 0 && length(atL1) > 0){
#        colz <- as.numeric(probe %in% bigD13) + 1
#        if(colz !=1) points(x = c(chrXaxt-1, chrXaxt+1), y = c(pos_c1[atC1], pos_l1[atL1]), lwd=colz, t='l', lty = c("dotted", "solid")[colz], col=c("purple", "gold")[colz])
#      }
    }
    points(x = rep(chrXaxt-1, length(pos_c1)), y = pos_c1, pch="-", col=c("gray", "white")[as.numeric(locs[on_c1,"chir1_Strand"])], cex=1)
    points(x = rep(chrXaxt-0, length(pos_c2)), y = pos_c2, pch="-", col=c("gray", "white")[as.numeric(locs[on_c2,"chir2_Strand"])], cex=1)
    points(x = rep(chrXaxt+1, length(pos_l1)), y = pos_l1, pch="-", col=c("gray", "white")[as.numeric(locs[on_l1,"lwlt01_Strand"])], cex=1)
    chrXaxt <- chrXaxt + 4
  }

  axis(1, at = seq(2, (length(chrs) * 4) + 1, 4), chrs)
  axis(2, at = seq(0, maxbp, 25000000), round(seq(0, maxbp, 25000000) / 1000000, 0), las=2)
  
  points(x = c(chrXaxt-5, chrXaxt-5), y = c(130000000, 150000000), lwd=5, t='l')
  text(chrXaxt-5, 120000000, "CHIR 1.0", srt = 90, cex=0.8)
  points(x = c(chrXaxt-4, chrXaxt-4), y = c(130000000, 150000000), lwd=5, t='l')
  text(chrXaxt-4, 120000000, "CHIR 2.0", srt = 90, cex=0.8)
  points(x = c(chrXaxt-3, chrXaxt-3), y = c(130000000, 150000000), lwd=5, t='l')
  text(chrXaxt-3, 120000000, "LWLT01", srt = 90, cex=0.8)

  #dev.off()
#}


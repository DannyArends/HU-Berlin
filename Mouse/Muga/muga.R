# MugaAnalysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")
chromosomes  <- c(1:19, "X", "Y", "M")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
mlength      <- max(chrInfo[,"Length"])

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
load("MM_snps.Rdata")

locusxdnaheader <- unlist(strsplit(readLines("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", n=16)[16],","))
locusxdna <- read.csv("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", header=FALSE, skip = 22)

colnames(locusxdna) <- c("Label","plateWell","Date","oligoPoolId","bundleId", "status", "Type", "Nas", locusxdnaheader[4:length(locusxdnaheader)])

locusxdna <- t(locusxdna[which(locusxdna[,"Type"] == "calls"),])
colnames(locusxdna) <- locusxdna[1,]
locusxdna <- locusxdna[-c(1:8), ]                                                     # Remove some unneeded headers

locusxdna[which(locusxdna == "U")] <- NA                                              # Use the correct NA values

MM_snps[,"SNP_ID"] <- gsub(".", "-", MM_snps[,"SNP_ID"], fixed = TRUE)                # Fix the . to - error in the MM_SNP
sum(rownames(locusxdna) %in% MM_snps[,"SNP_ID"])                                      # Should be: 77808
MM_snps <- MM_snps[match(rownames(locusxdna), MM_snps[,"SNP_ID"]),]                   # Sort MM_snps to match with our file

locusxdna[1:10,1:10]                                                                  # Convince ourselves that they match
MM_snps[1:10,1:5]                                                                     # Convince ourselves that they match

#cat("", file="Analysis/ArraySequences.fasta")                                         # Create a FASTA file to blast to the MM10 genome
#for(x in 1:nrow(MM_snps)){
#  cat(">", MM_snps[x,"SNP_ID"],"\n", file="Analysis/ArraySequences.fasta", append=TRUE)
#  cat(gsub("\\[.+\\]", "N", MM_snps[x,"Sequence"],perl=TRUE),"\n", file="Analysis/ArraySequences.fasta", append=TRUE)
#}

# Blast, we can use the database created for the Mouse Diversity Array
# blastn -task blastn -query ArraySequences.fasta -db E:\Mouse\DNA\DiversityArray\Analysis\Mus_musculus.GRCm38.74.dna.db -perc_identity 98 -outfmt 6 -evalue 0.1 -max_target_seqs 5 -out ProbeLocationBLAST.txt

# NOTE: we use the reported positions as by JAX

nGenotypes   <- apply(locusxdna, 1, function(x){ length(unique(na.omit(x))) })        # Number of unique genotypes
segregating  <- which(nGenotypes >= 2)                                                # Only segregating genotypes
genotypes    <- locusxdna[segregating,]                                               # Create the genotypes

imap <- MM_snps[segregating, c("Chr","Mb_NCBI38","cM")]                               # Get an initial map
rownames(imap) <- MM_snps[segregating, "SNP_ID"]
imap <- imap[order(as.numeric(imap$Mb_NCBI38)),]                                      # Sort it

map <- NULL
for(chr in chromosomes){ map <- rbind(map, imap[imap$Chr == chr,]); }                 # Correct the chromosome ordering
map$Mb_NCBI38 <- map$Mb_NCBI38 * 1000000                                              # We use basepairs, not megabases

genotypes <- genotypes[rownames(map),]                                                # Not all markers have a valid chromosome position, so take only the genotypes that survived

### Visualize the locations of the 'good' markers
plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="MegaMuga markers", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1,lwd=2)
  cnt <<- cnt + 1
})

aa <- apply(map, 1,function(x){
  yloc <- match(as.character(x["Chr"]), chromosomes); 
  if(as.character(x["Chr"]) == "X"){ cat(as.numeric(x["Mb_NCBI38"]),"\n"); }
  xloc <- as.numeric(x["Mb_NCBI38"])
  points(x=xloc, y=yloc + 0.1, pch="|", cex=0.9)
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)

write.table(map, "Analysis/mapRAW.txt", sep="\t")
write.table(genotypes, "Analysis/genotypesRAW.txt", sep="\t", quote = FALSE)

# Further work in phasing: beagle.R and grandparents.R
# QTL analysis: qtl_analysis.R

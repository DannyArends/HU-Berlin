#
# Get the ratio's and test for difference
#

# After SNP calling
setwd("/home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/Analysis")
matB6N    <- c("5070","5071","5072")
matBFMI   <- c("5073","5074","5075")
BFMI      <- c("4868", "5067")
B6N       <- c("5068", "5069")

samtools.exec       <- "samtools"

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

vcfdata         <- read.table("population.vcf", header = TRUE, colClasses="character")
write.table(cbind(vcfdata[,"CHROM"],vcfdata[,"POS"]),"SNPlocations.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

bamfiles        <- list.files(pattern = ".aligned.sorted.realigned.recalibrated.bam$")
reference       <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.dna"
reference.fa    <- paste0(reference, ".fa")        #!

for(x in bamfiles){
  # index the bamfile
  indexed.bamfile <-  gsub(".bam", ".bay", x)
  execute(paste0(samtools.exec, " index ", x), indexed.bamfile)
  cat("Mpileup of file:", x, "\n")
  execute(paste0("/home/neubert/Keane/samtools-1.2/samtools mpileup -uv -t DV -t DP -l SNPlocations.txt -f ", reference.fa, " ", x, " | bcftools call -c - > ",paste0(x,".vcf")),paste0(x,".vcf"))
}

# Load in the created VCF files, and extract the DP4 for all locations which exist across all files
vcffiles        <- list.files(pattern = ".aligned.sorted.realigned.bam.vcf$")
dp4datafull <- NA
for(x in vcffiles){
  cat(x," - ")
  vcfdata <- read.table(x, sep="\t",colClasses="character")
  vcfdata <- vcfdata[-which(grepl("INDEL", vcfdata[,8])),]
  write.table(cbind(vcfdata[, 1], vcfdata[, 2], vcfdata[, 4], "\t", gsub(",","\t", unlist(lapply(vcfdata[, 8], function(x){ sub(".*?DP4=(.*?);.*", "\\1", x)} )))), paste0(x,".dp4"), row.names=FALSE, col.names=FALSE, quote=FALSE)

  dp4data <- read.table(paste0(x,".dp4"),sep="\t")
  dp4data <- cbind(dp4data, dp4data[,2] + dp4data[,3])
  colnames(dp4data)[6] <- paste0(substr(x,1,4),"_Ref")
  dp4data <- cbind(dp4data, dp4data[,4] + dp4data[,5])
  colnames(dp4data)[7] <- paste0(substr(x,1,4),"_Alt")
  dp4data <- cbind(dp4data, Total = dp4data[,paste0(substr(x,1,4),"_Ref")] + dp4data[, paste0(substr(x,1,4),"_Alt")])
  if(is.na(dp4datafull)){
    dp4datafull <- dp4data[,c(1,6,7)]
  }else{
    dp4datafull <- dp4datafull[which(dp4datafull[,1] %in% dp4data[,1]),]
    dp4data <- dp4data[which(dp4data[,1] %in% dp4datafull[,1]),]
    dp4datafull <- cbind(dp4datafull, dp4data[,c(6,7)])
  }
  cat(nrow(dp4datafull),"\n")
}
write.table(dp4datafull, "allsamples.dp4",sep="\t",quote=FALSE, row.names=FALSE)

# Transfer to local machine

matB6N    <- c("5070","5071","5072"); matBFMI   <- c("5073","5074","5075") ; BFMI      <- c("4868", "5067") ; B6N       <- c("5068", "5069")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
reads <- read.table("allsamples.recal.dp4", sep="\t", header=TRUE, check.names=FALSE,row.names=1)       # load in the DP4 read count
reads <- reads[which(apply(reads,1,sum) >= 100),]                                                 # Filter for 100 reads across all samples

probs <- matrix(NA,nrow(reads), length(c(matB6N,matBFMI)),dimnames=list(rownames(reads), c(matB6N,matBFMI)))      # Raw ChiSquare probabilities
ratios <- matrix(NA,nrow(reads), length(c(matB6N,matBFMI)),dimnames=list(rownames(reads), c(matB6N,matBFMI)))     # ref / total ratios
chisq <- matrix(NA,nrow(reads), length(c(matB6N,matBFMI)),dimnames=list(rownames(reads), c(matB6N,matBFMI)))      # ChiSquare test statistic
for(sname in c(matB6N,matBFMI)){
  for(x in rownames(reads)){
    totalreads <- as.numeric(reads[x, paste0(sname,"_Ref")]) + as.numeric(reads[x, paste0(sname,"_Alt")])
    observed <- c(as.numeric(reads[x, paste0(sname,"_Ref")]), as.numeric(reads[x, paste0(sname,"_Alt")]))
    if(totalreads > 5){
      ratios[x, sname] <- round(as.numeric(reads[x, paste0(sname,"_Alt")]) / totalreads,3) * 100
      probs[x, sname] <- chisq.test(observed)$p.value
      chisq[x, sname] <- sign(ratios[x, sname] - 50) * sqrt(as.numeric(unlist(chisq.test(observed))["statistic.X-squared"]))
    }
  }
  cat("Done", sname,"\n")
}
write.table(ratios, "F1_ChiSquare.ratios.txt", sep="\t")
write.table(probs,  "F1_ChiSquare.probs.txt", sep="\t")
write.table(chisq,  "F1_ChiSquare.adjusted.txt", sep="\t")

# ChiSquare test of parental allele (Ref, Hetro, Alt)
pAllele <- matrix(NA,nrow(reads), length(c(B6N,BFMI)), dimnames=list(rownames(reads), c(B6N,BFMI)))
for(sname in c(B6N, BFMI)){
  for(x in rownames(reads)){
    totalreads <- as.numeric(reads[x, paste0(sname,"_Ref")]) + as.numeric(reads[x, paste0(sname,"_Alt")])
    if(totalreads > 5){
      observed <- c(as.numeric(reads[x, paste0(sname,"_Ref")]), as.numeric(reads[x, paste0(sname,"_Alt")]))
      i <- which.max(c(chisq.test(rbind(observed, c(totalreads, 1)))$p.value, chisq.test(observed)$p.value, chisq.test(rbind(observed, c(1, totalreads)))$p.value))
      pAllele[x, sname] <- c("Ref", "Hetro", "Alt")[i]
    }
  }
  cat("Done", sname,"\n")
}
write.table(pAllele, "parentalAlleles.txt", sep="\t")

# Load in the GTF data, and make a lookup table for all exons in the mouse genome
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
GTF <- read.table("GTF/Mus_musculus.GRCm38.81.gtf", sep="\t")                                                       # Gene models
EXONS <- GTF[which(GTF[,3]=="exon"),]

Gene_Exon <- t(apply(EXONS, 1, function(x){
  c(sub(".*?gene_id (.*?);.*", "\\1", x[9]), sub(".*?exon_id (.*?);.*", "\\1", x[9]), sub(".*?transcript_id (.*?);.*", "\\1", x[9]), sub(".*?exon_number (.*?);.*", "\\1", x[9]), sub(".*?gene_name (.*?);.*", "\\1", x[9]))
}))
Gene_Exon <- cbind(Gene_Exon, EXONS[, c(1, 4, 5)])
colnames(Gene_Exon) <- c("gene_id", "exon_id", "transcript_id", "exon_number", "gene_name", "chr", "start", "end")

# Load in the created files
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
probs   <- read.table("F1_ChiSquare.probs.txt", sep="\t",row.names=1,header=TRUE, check.names=FALSE)
ratios  <- read.table("F1_ChiSquare.ratios.txt", sep="\t",row.names=1,header=TRUE, check.names=FALSE)
chisq   <- read.table("F1_ChiSquare.adjusted.txt", sep="\t",row.names=1,header=TRUE, check.names=FALSE)
pAllele <- read.table("parentalAlleles.txt", sep="\t",row.names=1,header=TRUE, check.names=FALSE,colClasses="character")

# Annotate the SNPs using GTF information
results <- NULL
for(x in rownames(pAllele)){
  snpinfo  <- strsplit(x, " ")[[1]]
  snp.chr  <-  snpinfo[1]
  snp.loc  <-  as.numeric(snpinfo[2])
  cat("Going to look for a gene at:", snp.chr, snp.loc)
  idx <- which(Gene_Exon[,"chr"] == snp.chr & Gene_Exon[,"start"] <= snp.loc & Gene_Exon[,"end"] >= snp.loc)
  cat("", length(idx),"")
  if(length(idx) == 0){
    mrow <- c(x, rep("", ncol(Gene_Exon)), pAllele[x,], ratios[x,], chisq[x,], probs[x,])
    cat(length(mrow))
    results <- rbind(results, mrow)
  }else if(length(idx) == 1){
    mrow <- c(x, unlist(lapply(Gene_Exon[idx,],as.character)), pAllele[x,], ratios[x,], chisq[x,], probs[x,])
    cat(length(mrow))
    results <- rbind(results, mrow)
  }else{
    mrow <- cbind(x, apply(Gene_Exon[idx,],2,as.character))
    pRow <- matrix(unlist(rep(pAllele[x,], nrow(mrow))), nrow(mrow), ncol(pAllele),byrow=TRUE, dimnames = list(1:nrow(mrow),colnames(pAllele[x,])))
    rRow <- matrix(unlist(rep(ratios[x,],  nrow(mrow))), nrow(mrow), ncol(ratios), byrow=TRUE, dimnames = list(1:nrow(mrow),colnames(ratios[x,])))
    cRow <- matrix(unlist(rep(chisq[x,],   nrow(mrow))), nrow(mrow), ncol(chisq),  byrow=TRUE, dimnames = list(1:nrow(mrow),colnames(chisq[x,])))
    oRow <- matrix(unlist(rep(probs[x,],   nrow(mrow))), nrow(mrow), ncol(probs),  byrow=TRUE, dimnames = list(1:nrow(mrow),colnames(probs[x,])))
    mrow <- cbind(mrow, pRow, rRow, cRow, oRow)
    cat(ncol(mrow))
    results <- rbind(results, mrow)
  }
  cat("\n")
}
rownames(results) <- 1:nrow(results)
write.table(results,"annotatedSNPsChiSq.txt", sep = "\t")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/GTF/")
canonical <- read.table("mouse_canonical_transcripts_v81.fa",sep="\t", check.names=FALSE,colClasses="character")
canonical <- rbind(canonical, c("",""))                                                         # Empty transcript for genes not in exons

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
results <- read.table("annotatedSNPsChiSq.txt",sep="\t", header=TRUE, check.names=FALSE)
results <- results[which(results[,"transcript_id"] %in% canonical[,2]),]                        # Only take the canonical transcripts
colnames(results)[20:25] <- paste0(colnames(results)[20:25],"_ChiSq")                           # Add the correct headers
colnames(results)[26:31] <- paste0(colnames(results)[26:31],"_P")                               # Add the correct headers


isEqual <- function(x){
  if(is.na(x[1]) || is.na(x[2])) return(TRUE)
  return(x[1] == x[2])
}

# Filer out where the 2 BFMIs or the 2 B6N individuals are different
equalBFMI <- apply(results[,BFMI], 1, isEqual)
results <- results[equalBFMI, ]

equalB6N <- apply(results[,B6N], 1, isEqual)
results <- results[equalB6N, ]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After basic QC there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")
#After basic QC there are 45300 SNPs, there are 16099 not in genes, and 29201 SNPs in 5265 genes

BFMIalleles <- apply(results[,BFMI],1,function(x){
  if(is.na(x[1])) return(x[2])
  return(x[1])
})

B6Nalleles <- apply(results[,B6N],1,function(x){
  if(is.na(x[1])) return(x[2])
  return(x[1])
})

onAutosome <- which(results[,"chr"] != "X" & results[,"chr"] != "Y" & results[,"chr"] != "MT")
results <- results[onAutosome,]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting autosomes there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")
#After selecting autosomes there are 44024 SNPs, there are 16099 not in genes, and 27925 SNPs in 5082 genes

# At which SNPs do the father and the mother differ ?
equalto <- rep(FALSE, nrow(results))
for(x in 1:nrow(results)){
  equalto[x] <- (BFMIalleles[x] == B6Nalleles[x])
}

results <- results[which(!equalto),]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting BFMI != B6N there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")
##After selecting BFMI != B6N there are 25515 SNPs, there are 6759 not in genes, and 18756 SNPs in 3707 genes

plot(apply(results[,paste0(matB6N,"_ChiSq")],1,sum) - apply(results[,paste0(matBFMI,"_ChiSq")],1,sum))

matB6N.sum <- apply(results[,paste0(matB6N,"_ChiSq")],1,sum,na.rm=TRUE)
matBFMI.sum <- apply(results[,paste0(matBFMI,"_ChiSq")],1,sum,na.rm=TRUE)


oldResults <- read.table("Table S2.txt",sep="\t", header = TRUE)
snpIDnames <- paste(oldResults[,"CHROM"], oldResults[,"POS"], oldResults[,"REF"],"")
snpShared <- snpIDnames[which(snpIDnames %in% as.character(results[,1]))]
cat("SNPs shared between 100 % analysis and chiSq =", round(length(snpShared) / length(snpIDnames),3) * 100,"%\n")

write.table(results[which(as.character(results[,1]) %in% snpShared),],"ChiSquareTableS2.txt",sep="\t",row.names=FALSE,quote=FALSE)

#
# Find things which are different, use permutation on our data
#

randomScoresSum <- NULL
randomScores <- NULL
for(x in 1:10000){
  i <- sample(nrow(results), 1)
  randomScoresSum <- c(randomScoresSum, max(abs(c(matB6N.sum[i],matBFMI.sum[i])),na.rm=TRUE))
  randomScores <- c(randomScores, max(abs(results[i,c(paste0(matB6N,"_ChiSq"), paste0(matBFMI,"_ChiSq"))]),na.rm=TRUE))
}
threshold.a5S <- sort(randomScoresSum)[length(randomScoresSum) * .95] ; threshold.a5 <- sort(randomScores)[length(randomScores) * .95]
threshold.a1S <- sort(randomScoresSum)[length(randomScoresSum) * .99] ; threshold.a1 <- sort(randomScores)[length(randomScores) * .99]

showsASE <- results[unique(c(which(matB6N.sum > threshold.a5S | matB6N.sum < -threshold.a5S), which(matBFMI.sum > threshold.a5S | matBFMI.sum < -threshold.a5S))),]
write.table(showsASE, "ChiSquareAllSupplTable.txt",sep="\t",row.names=FALSE,quote=FALSE)

diffs <- NULL
ratios <- NULL
for(x in 1:nrow(results)){
  mSq.b6n  <- mean(as.numeric(results[x, paste0(matB6N,"_ChiSq")]),na.rm=TRUE)
  mRa.b6n  <- mean(as.numeric(results[x, paste0(matB6N)]),na.rm=TRUE)
  mSq.bfmi <- mean(as.numeric(results[x, paste0(matBFMI,"_ChiSq")]),na.rm=TRUE)
  mRa.bfmi <- mean(as.numeric(results[x, paste0(matBFMI)]),na.rm=TRUE)
  mdiff <- abs(mSq.b6n - mSq.bfmi)
  mrat  <- abs(mRa.b6n - mRa.bfmi)
  cat(x,  mRa.b6n, mRa.bfmi, mrat, mdiff,"\n")
  diffs <- c(diffs,mdiff)
  ratios <- c(ratios, mrat)
}

getAdjustedChiSq <- function(ratio, chiSq){
  if(is.na(ratio)) return(NA)
  if(is.na(chiSq)) return(NA)
  if(ratio < 50) return(-chiSq)
  return(chiSq)
}

adjmBFMI <- NULL
adjmB6N<- NULL
for(x in 1:nrow(results)){
  mBFMI1.adjChiSq <- getAdjustedChiSq(results[x, matBFMI[1]], results[x, paste0(matBFMI[1],"_ChiSq")])
  mBFMI2.adjChiSq <- getAdjustedChiSq(results[x, matBFMI[2]], results[x, paste0(matBFMI[2],"_ChiSq")])
  mBFMI3.adjChiSq <- getAdjustedChiSq(results[x, matBFMI[3]], results[x, paste0(matBFMI[3],"_ChiSq")])
  adjmBFMI <- c(adjmBFMI, sum(c(mBFMI1.adjChiSq, mBFMI2.adjChiSq, mBFMI3.adjChiSq),na.rm=TRUE))

  mB6N1.adjChiSq <- getAdjustedChiSq(results[x, matB6N[1]], results[x, paste0(matB6N[1],"_ChiSq")])
  mB6N2.adjChiSq <- getAdjustedChiSq(results[x, matB6N[2]], results[x, paste0(matB6N[2],"_ChiSq")])
  mB6N3.adjChiSq <- getAdjustedChiSq(results[x, matB6N[3]], results[x, paste0(matB6N[3],"_ChiSq")])
  adjmB6N <- c(adjmB6N, sum(c(mB6N1.adjChiSq, mB6N2.adjChiSq, mB6N3.adjChiSq),na.rm=TRUE))
}


good <- NULL

for(x in 1:nrow(results)){
  ase.bfmi <- all(results[x, paste0(matBFMI,"_P")] < 0.05,na.rm=TRUE)
  ase.b6n <- all(results[x, paste0(matB6N,"_P")] < 0.05,na.rm=TRUE)
  if(sum(c(ase.bfmi, ase.b6n)) == 1){ # One side shows ASE
    cat(x, ase.bfmi, ase.b6n,"\n")
    good <- c(good, x)
  }
  if(sum(c(ase.bfmi, ase.b6n)) == 2){
    # Both sides show ASE
    b6.d <- all(results[x, matB6N] <= 50,na.rm=TRUE) ; b6.u <- all(results[x, matB6N] >= 50,na.rm=TRUE)
    bf.d <- all(results[x, matBFMI] <= 50,na.rm=TRUE) ; bf.u <- all(results[x, matBFMI] >= 50,na.rm=TRUE)
    cat(x, b6.d, b6.u, bf.d, bf.u,"\n")
  }
}

results2 <- results[which(abs(adjmBFMI-adjmB6N) > 20),]


chrs <- unlist(lapply(strsplit(as.character(results[,1])," "),"[",1))
onAuto <- which(chrs != "X" & chrs != "MT")

hist(abs(adjmBFMI-adjmB6N)[onAuto])

results2 <- results[which(ratios > 40),]







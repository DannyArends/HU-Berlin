#
# Quality control (canonical transcripts), after whcih we do filtering for candidates (B66N != BFMI)
# After which we get the alt/(ref+alt) ratio's and test for difference in ChiSquare values using permutations
#

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
results <- results[-which(apply(results[,c(BFMI,B6N)],1,function(x){return(any(x == "Hetro"))})), ] #Throw away SNPs at which the parental strains are hetrozygous


snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting BFMI != B6N there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")
#After selecting BFMI != B6N there are 24944 SNPs, there are 6565 not in genes, and 18379 SNPs in 3630 genes

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
randomScoresDiff <- NULL
randomScores <- NULL
for(x in 1:10000){
  i <- sample(nrow(results), 1)
  randomScoresSum <- c(randomScoresSum, max(abs(c(matB6N.sum[i],matBFMI.sum[i])),na.rm=TRUE))
  randomScoresDiff <- c(randomScoresSum, max(abs(matB6N.sum[i] - matBFMI.sum[i]),na.rm=TRUE))
  randomScores <- c(randomScores, max(abs(results[i,c(paste0(matB6N,"_ChiSq"), paste0(matBFMI,"_ChiSq"))]),na.rm=TRUE))
}
threshold.a5S <- sort(randomScoresSum)[length(randomScoresSum) * .95] ; threshold.a5 <- sort(randomScores)[length(randomScores) * .95] ; threshold.a5d <- sort(randomScoresDiff)[length(randomScoresDiff) * .95]
threshold.a1S <- sort(randomScoresSum)[length(randomScoresSum) * .99] ; threshold.a1 <- sort(randomScores)[length(randomScores) * .99] ; threshold.a1d <- sort(randomScoresDiff)[length(randomScoresDiff) * .99]

results <- results[unique(c(which(matB6N.sum > threshold.a5S | matB6N.sum < -threshold.a5S), which(matBFMI.sum > threshold.a5S | matBFMI.sum < -threshold.a5S), which(matB6N.sum-matBFMI.sum > threshold.a5d | matB6N.sum-matBFMI.sum < -threshold.a5d))),]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting for ASE there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")

### Add the mgi descriptions to the genes.
ensembleIDs <- as.character(unique(showsASE[,"gene_id"]))
library(biomaRt)
mart      <- useMart("ensembl", "mmusculus_gene_ensembl")
descriptions  <- getBM(values = ensembleIDs, filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_description"), mart = mart)

results <- cbind(results, mgi_description = NA)
for(x in 1:nrow(results)){
  ix <- which(descriptions[,"ensembl_gene_id"]  == results[x,"gene_id"])
  if(length(ix) == 1)results[x, "mgi_description"] <- descriptions[ix, "mgi_description"]
}

results <- cbind(results, diff = NA)
results <- cbind(results, absdiff = NA)
for(x in 1:nrow(results)){
  results[x,"diff"] <- sum(results[x,paste0(matB6N,"_ChiSq")],na.rm=TRUE) - sum(results[x,paste0(matBFMI,"_ChiSq")],na.rm=TRUE)
  results[x,"absdiff"] <- abs(results[x,"diff"])
}

results <- cbind(results, meanRatioMatB6N = apply(results[,matB6N],1,mean,na.rm=TRUE))
results <- cbind(results, sdRatioMatB6N = apply(results[,matB6N],1,sd,na.rm=TRUE))
results <- cbind(results, meanRatioMatBFMI = apply(results[,matBFMI],1,mean,na.rm=TRUE))
results <- cbind(results, sdRatioMatBFMI = apply(results[,matBFMI],1,sd,na.rm=TRUE))
results <- cbind(results, diffRatio = abs(results[,"meanRatioMatB6N"] - results[,"meanRatioMatBFMI"]))

## High standard deviations  > 20 in one of the groups and you are out
results <- results[-which(results[,"sdRatioMatB6N"] > 20 | results[,"sdRatioMatBFMI"] > 20),]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting for ASE there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")
write.table(results, "ChiSquareAllSupplTable.txt",sep="\t",row.names=FALSE,quote=FALSE)

# Filter results:
# - False positives (Might be interesting because of the ifference in expression reflected in the ChiSq scores
# - True positives
# - False Negatives / Borderline

write.table(results, "ChiSquareFilteredSupplTable.txt",sep="\t",row.names=FALSE,quote=FALSE)


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







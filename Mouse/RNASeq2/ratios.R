#
# Quality control (canonical transcripts), after whcih we do filtering for candidates (B66N != BFMI)
# After which we get the alt/(ref+alt) ratio's and test for difference in ChiSquare values using permutations
#

library(biomaRt)

matB6N    <- c("5070","5071","5072"); matBFMI   <- c("5073","5074","5075") ; BFMI      <- c("4868", "5067") ; B6N       <- c("5068", "5069")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
canonical <- read.table("GTF/mouse_canonical_transcripts_v81.fa",sep="\t", check.names=FALSE,colClasses="character")    # List of canonical transcripts
canonical <- rbind(canonical, c("",""))                                                                                 # Empty transcript for genes not in exons

results <- read.table("ReAnalysisSNPs/annotatedSNPsChiSq.txt",sep="\t", header=TRUE, check.names=FALSE)

colnames(results)[20:25] <- paste0(colnames(results)[20:25],"_ChiSq")                                                   # Add the correct headers
colnames(results)[26:31] <- paste0(colnames(results)[26:31],"_P")                                                       # Add the correct headers

singletons <- which(!results[,1] %in% results[duplicated(results[,1]),1])                                               # SNPs only in the results once

passed <- results[singletons,]                                                          # If SNPs occur only once we pass the SNP on to the next round of analysis
results <- results[-singletons,]                                                        # SNPs that occur multiple times, we need to check which we need to take (canonical, etc)

for(x in 1:nrow(results)){
  snpID <- results[x, 1]
  if(!(snpID %in% passed[,1])){                                                         # If this SNP is not yet passed
    ix <- which(results[,1] == snpID)
    if(length(ix) == 1){                                                                # This SNP only occurs once (Should not occur since we passed all single SNPs beforehand
      passed <- rbind(passed, results[x,])
    }else{                                                                              # SNPs is in multiple exons / genes
      msubset <- results[ix,]
      genes <- unique(msubset[, "gene_id"])                                             # Get all the genes this SNP is in
      for(gene in genes){
        iy <- which(canonical[,1] == gene)
        if(length(iy) == 1){
          canonicalTransciptName <- canonical[iy,2]
          Tsubset <- msubset[which(msubset[,"gene_id"] == gene),]
          iz <- which(Tsubset[,"transcript_id"] == canonicalTransciptName)
          if(length(iz) == 1){                                                          # SNP is in an exon of the canonical transcript
            passed <- rbind(passed, Tsubset[iz,])
          }else if(length(iz) == 0){                                                    # SNP is not in an exon of the canonical transcript
            passed <- rbind(passed, Tsubset[1,])
          }else{                                                                        # SNP is in 2 exons, which are both canonical transcript for this gene ?! 
            stop("Should not occur, canonical transcript 2 times")                      # This should not occur
          }
        }else if(length(iy)  == 0){
          Tsubset <- msubset[which(msubset[,"gene_id"] == gene),]
          passed <- rbind(passed, Tsubset[1,])
        }else{
          stop("Should not occur, duplicate gene in canonical list")
        }
      }
    }
  }
  cat(x,",", nrow(results),"passed:",nrow(passed),"\n")
}
results <- passed                                                                           # SNPs that passed are now the results
write.table(results, "ReAnalysisSNPs/annotatedSNPsChiSq_filtered.txt", sep = "\t")          # Save the annotated and duplicate filtered results

# Reload the annotated SNPs form here
library(biomaRt)
matB6N    <- c("5070","5071","5072"); matBFMI   <- c("5073","5074","5075") ; BFMI      <- c("4868", "5067") ; B6N       <- c("5068", "5069")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
results <- read.table("ReAnalysisSNPs/annotatedSNPsChiSq_filtered.txt", sep = "\t", stringsAsFactors =FALSE,check.names=FALSE)
colnames(results)[20:25] <- paste0(colnames(results)[20:25],"_ChiSq")                       # Add the correct headers
colnames(results)[26:31] <- paste0(colnames(results)[26:31],"_P")                           # Add the correct headers

# Do some additional annotation of SNPs in introns of genes
mart      <- useMart("ensembl", "mmusculus_gene_ensembl")
allgenes  <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), mart = mart)

for(x in 1:nrow(results)){
  if(results[x, "gene_id"] == ""){
    msplit <- strsplit(as.character(results[x, 1])," ")[[1]]
    cat(x, msplit)
    ix <- which(allgenes[,"chromosome_name"] == msplit[1] & allgenes[,"start_position"] <= as.numeric(msplit[2]) & allgenes[,"end_position"] >= as.numeric(msplit[2]))
    if(length(ix) == 1){
      #cat(", in gene:", allgenes[ix,"ensembl_gene_id"], ":", allgenes[ix,"start_position"], "<", as.numeric(msplit[2]) , "<", allgenes[ix,"end_position"])
      results[x,"gene_id"]   <- as.character(allgenes[ix,"ensembl_gene_id"])
      results[x,"chr"]       <- as.character(allgenes[ix,"chromosome_name"])
      results[x,"start"]     <- as.character(allgenes[ix,"start_position"])
      results[x,"end"]       <- as.character(allgenes[ix,"end_position"])
    }
    cat("\n")
  }
}
write.table(results, "ReAnalysisSNPs/annotatedSNPsChiSq_filtered_inGene.txt", sep = "\t")

# Reload the annotated SNPs form here
matB6N    <- c("5070","5071","5072"); matBFMI   <- c("5073","5074","5075") ; BFMI      <- c("4868", "5067") ; B6N       <- c("5068", "5069")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
results <- read.table("ReAnalysisSNPs/annotatedSNPsChiSq_filtered_inGene.txt", sep = "\t",header=TRUE, check.names=FALSE,stringsAsFactors=FALSE)

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
#After basic QC there are 49120 SNPs, there are 7622 not in genes, and 41498 SNPs in 6479 genes

chrs <- unlist(lapply(lapply(results[,1],strsplit," "),function(x){return(x[[1]][1])}))

onAutosome <- which(chrs != "X" & chrs != "Y" & chrs != "MT")
results <- results[onAutosome,]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting autosomes there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")
#After selecting autosomes there are 44024 SNPs, there are 16099 not in genes, and 27925 SNPs in 5082 genes

BFMIalleles <- apply(results[,BFMI],1,function(x){ if(is.na(x[1])) return(x[2]) ; return(x[1]) })
B6Nalleles <- apply(results[,B6N],1,function(x){ if(is.na(x[1])) return(x[2]) ; return(x[1]) })

# At which SNPs do the father and the mother differ ?
equalto <- rep(FALSE, nrow(results))
for(x in 1:nrow(results)){
  equalto[x] <- (BFMIalleles[x] == B6Nalleles[x])
}

results <- results[which(!equalto),]
results <- results[-which(apply(results[,c(BFMI,B6N)],1,function(x){return(any(x == "Hetro"))})), ] #Throw away SNPs at which a parental strain is heterozygous

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting BFMI != B6N there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")
#After selecting BFMI != B6N there are 24944 SNPs, there are 6565 not in genes, and 18379 SNPs in 3630 genes

results[which(results[,1] == "1 9553659 G "),1:15]

matB6N.sum <- apply(results[,paste0(matB6N,"_ChiSq")],1,sum,na.rm=TRUE)
matBFMI.sum <- apply(results[,paste0(matBFMI,"_ChiSq")],1,sum,na.rm=TRUE)

plot(matB6N.sum, main="Summed Chi2 scores matB6N group")
plot(matBFMI.sum, main="Summed Chi2 scores matBFMI group")
plot(matBFMI.sum - matB6N.sum, main="Difference in Chi2 scores matBFMI - matB6N")

# Compare back to our old results (which only looked for mono-allelic expression)
oldResults <- read.table("ReAnalysisSNPs/Table S2.txt",sep="\t", header = TRUE)
snpIDnames <- paste(oldResults[,"CHROM"], oldResults[,"POS"], oldResults[,"REF"],"")
snpShared <- snpIDnames[which(snpIDnames %in% as.character(results[,1]))]
cat("SNPs shared between 100 % analysis and chiSq =", round(length(snpShared) / length(snpIDnames),3) * 100,"%\n")

# Write the overlap to a separate file
write.table(results[which(as.character(results[,1]) %in% snpShared),],"ReAnalysisSNPs/ChiSquareTableS2.txt",sep="\t",row.names=FALSE,quote=FALSE)

## Permutation to find which delta Chi2 scores are significantly different
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

significantASE <- unique(c(which(matB6N.sum > threshold.a5S | matB6N.sum < -threshold.a5S),                                 # Significant ASE in matB6N 
                           which(matBFMI.sum > threshold.a5S | matBFMI.sum < -threshold.a5S),                               # Significant ASE in matBFMI 
                           which(matB6N.sum-matBFMI.sum > threshold.a5d | matB6N.sum-matBFMI.sum < -threshold.a5d)))        # Significant difference in ASE
results <- results[significantASE,]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting for ASE there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")

### Add the mgi descriptions to the genes.
ensembleIDs <- as.character(unique(results[,"gene_id"]))
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

## If we find a high standard deviations of the alt/(alt+ref) ratio ( > 20 ) in one of the groups the SNP is not considered
results <- results[-which(results[,"sdRatioMatB6N"] > 20 | results[,"sdRatioMatBFMI"] > 20),]

snp.unique    <- length(unique(results[,1]))
snp.outside   <- nrow(results[which(results[,2] == ""),])
snp.ingenes   <- length(unique(results[which(results[,2] != ""),1]))
genes.unique  <- length(unique(results[which(results[,2] != ""),2]))

cat("After selecting for ASE and filtering for concordance there are", snp.unique, "SNPs, there are", snp.outside, "not in genes, and", snp.ingenes, "SNPs in", genes.unique ,"genes\n")

# Add expression data from: BFMI_RPKM_Qnorm_ANN_AddDom_plusLog2.txt
expfile <- read.table("Analysis/BFMI_RPKM_Qnorm_ANN_AddDom_plusLog2.txt", sep="\t",header=TRUE,stringsAsFactors=FALSE)
results <- cbind(results, Mean.BFMI860.12xB6N.L_log = NA)
results <- cbind(results, Mean.B6NxBFMI860.12.L_log = NA)
results <- cbind(results, tTest_F1_log = NA)
results <- cbind(results, Ratio_F1_log = NA)

for(x in 1:nrow(results)){
  l <- which(rownames(expfile) == results[x,"gene_id"]);
  if(length(l) == 1){
    results[x,"Mean.BFMI860.12xB6N.L_log"] <- expfile[l,"Mean.BFMI860.12xB6N.L_log"]
    results[x,"Mean.B6NxBFMI860.12.L_log"] <- expfile[l,"Mean.B6NxBFMI860.12.L_log"]
    results[x,"tTest_F1_log"] <- expfile[l,"tTest_F1_log"]
    results[x,"Ratio_F1_log"] <- expfile[l,"Ratio_F1_log"]
  }
}

splitup <- lapply(as.character(results[,1]), strsplit, " ")
chrs <- unlist(lapply(splitup,function(x){return(x[[1]][1])}))
pos <- unlist(lapply(splitup,function(x){return(x[[1]][2])}))
ref <- unlist(lapply(splitup,function(x){return(x[[1]][3])}))

results[,"Chr"] <- chrs
results <- cbind(results, Pos = pos)
results <- cbind(results, Ref = ref)

populationVCF <- read.table("ReAnalysisSNPs/population.vcf",sep="\t",header=TRUE)
idsInPOP <- paste0(populationVCF[,1]," ", populationVCF[,2]," ", populationVCF[,4]," ", sep="")

results <- cbind(results, Alt = NA)
for(x in 1:nrow(results)){
  ii <- which(idsInPOP == results[x,1])
  if(length(ii) == 1) results[x,"Alt"] <- as.character(populationVCF[ii,"ALT"])
}


# TODO: Add the generegion (intron, exon, intergenic)


write.table(results, "ReAnalysisSNPs/ChiSquareAllSupplTable.txt",sep="\t",row.names=FALSE,quote=FALSE)

plot(x = results[,"absdiff"], y = results[,"diffRatio"],xlab="∑ abs(Ref/(Ref+alt) - 0.5) * sqrt(chi2)", ylab = "Difference in ratio Ref/(Ref+alt) matBFMI vs matB6N",pch=19, main="ASE between matBFMI and matB6N")
abline(h=20)
abline(v=threshold.a5d)

# Filter results:
# - False positives (Might be interesting because of the ifference in expression reflected in the ChiSq scores
# - True positives
# - False Negatives / Borderline
write.table(results[which(results[,"absdiff"] > threshold.a5d & results[,"diffRatio"] > 20),], "ReAnalysisSNPs/ChiSquareFilteredSupplTable_ratioAbove20_Significant0.05.txt",sep="\t",row.names=FALSE,quote=FALSE)
write.table(results[which(results[,"diffRatio"] > 20),], "ReAnalysisSNPs/ChiSquareFilteredSupplTable_ratioAbove20.txt",sep="\t",row.names=FALSE,quote=FALSE)

#
# Chi Square calculation and annotation of the results using the GTF file
#

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

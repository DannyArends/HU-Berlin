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

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
reads <- read.table("allsamples.recal.dp4", sep="\t", header=TRUE, check.names=FALSE,row.names=1)       # load in the DP4 read count
reads <- reads[which(apply(reads,1,sum) >= 100),]                                                 # Filter for 100 reads across all samples

# ChiSquare probabilities
probs <- matrix(NA,nrow(reads), length(c(matB6N,matBFMI)),dimnames=list(rownames(reads), c(matB6N,matBFMI)))
for(sname in c(matB6N,matBFMI)){
  for(x in rownames(reads)){
    totalreads <- as.numeric(reads[x, paste0(sname,"_Ref")]) + as.numeric(reads[x, paste0(sname,"_Alt")])
    observed <- c(as.numeric(reads[x, paste0(sname,"_Ref")]), as.numeric(reads[x, paste0(sname,"_Alt")]))
    if(totalreads > 5) probs[x, sname] <- chisq.test(observed)$p.value
  }
  cat("Done", sname,"\n")
}
write.table(probs, "F1_ChiSquare.probs.txt", sep="\t")

# ref / total ratios
ratios <- matrix(NA,nrow(reads), length(c(matB6N,matBFMI)),dimnames=list(rownames(reads), c(matB6N,matBFMI)))
for(sname in c(matB6N,matBFMI)){
  for(x in rownames(reads)){
    totalreads <- as.numeric(reads[x, paste0(sname,"_Ref")]) + as.numeric(reads[x, paste0(sname,"_Alt")])
    if(totalreads > 5) ratios[x, sname] <- as.numeric(reads[x, paste0(sname,"_Alt")]) / totalreads
  }
  cat("Done", sname,"\n")
}
write.table(ratios, "F1_ChiSquare.ratios.txt", sep="\t")

# ChiSquare test statistic
chisq <- matrix(NA,nrow(reads), length(c(matB6N,matBFMI)),dimnames=list(rownames(reads), c(matB6N,matBFMI)))
for(sname in c(matB6N,matBFMI)){
  for(x in rownames(reads)){
    totalreads <- as.numeric(reads[x, paste0(sname,"_Ref")]) + as.numeric(reads[x, paste0(sname,"_Alt")])
    observed <- c(as.numeric(reads[x, paste0(sname,"_Ref")]), as.numeric(reads[x, paste0(sname,"_Alt")]))
    if(totalreads > 5) chisq[x, sname] <- unlist(chisq.test(observed))["statistic.X-squared"]
  }
  cat("Done", sname,"\n")
}
write.table(chisq, "F1_ChiSquare.txt", sep="\t")

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
    #    cat(x, reads[x,1], c("Ref", "Hetro", "Alt")[i], "\n")
    #  cat(x, reads[x,1], chisq.test(observed, p=c(.99, .01))$p.value, chisq.test(observed)$p.value, chisq.test(observed, p=c(.01, .99))$p.value,"\n")
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
chisq   <- read.table("F1_ChiSquare.txt", sep="\t",row.names=1,header=TRUE, check.names=FALSE)
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
    mrow <- c(x, rep("", ncol(Gene_Exon)), pAllele[x,], round(ratios[x,],3) * 100, chisq[x,], probs[x,])
    cat(length(mrow))
    results <- rbind(results, mrow)
  }else if(length(idx) == 1){
    mrow <- c(x, unlist(lapply(Gene_Exon[idx,],as.character)), pAllele[x,], round(ratios[x,],3) * 100, chisq[x,], probs[x,])
    cat(length(mrow))
    results <- rbind(results, mrow)
  }else{
    mrow <- cbind(x, apply(Gene_Exon[idx,],2,as.character))
    pRow <- matrix(unlist(rep(pAllele[x,], nrow(mrow))),nrow(mrow),ncol(pAllele),byrow=TRUE,dimnames = list(1:nrow(mrow),colnames(pAllele[x,])))
    rRow <- matrix(unlist(rep(round(ratios[x,],3) * 100, nrow(mrow))),nrow(mrow),ncol(ratios),byrow=TRUE,dimnames = list(1:nrow(mrow),colnames(ratios[x,])))
    cRow <- matrix(unlist(rep(chisq[x,], nrow(mrow))),nrow(mrow),ncol(chisq),byrow=TRUE,dimnames = list(1:nrow(mrow),colnames(chisq[x,])))
    oRow <- matrix(unlist(rep(probs[x,], nrow(mrow))),nrow(mrow),ncol(probs),byrow=TRUE,dimnames = list(1:nrow(mrow),colnames(probs[x,])))
    mrow <- cbind(mrow, pRow, rRow, cRow, oRow)
    cat(ncol(mrow))
    results <- rbind(results, mrow)
  }
  cat("\n")
}
rownames(results) <- 1:nrow(results)
write.table(results,"annotatedSNPsChiSq.txt", sep = "\t")

# At which SNPs do the father and the mother differ ?
matB6Ndiff <- NULL
matBFMIdiff <- NULL
for(x in rownames(pAllele)){
  if(!(is.na(pAllele[x,"4868"]) || is.na(pAllele[x,"5069"]))){
    if(pAllele[x,"4868"] != pAllele[x,"5069"]){ #cat(x, pAllele[x, c("5069", "4868")], probs[x, matB6N],"\n")
    matB6Ndiff <- c(matB6Ndiff, x)
    }
  }
  if(!(is.na(pAllele[x,"5068"]) || is.na(pAllele[x,"5067"]))){
    if(pAllele[x,"5068"] != pAllele[x,"5067"]){ #cat(x, pAllele[x, c("5068", "5067")], probs[x, matBFMI],"\n")
      matBFMIdiff <- c(matBFMIdiff, x)
    }
  }
}
useableProbes <- unique(c(matBFMIdiff, matB6Ndiff))




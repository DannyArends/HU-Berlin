# snpanalysis.R - Analyze the SNPs and indels called by the GenomeAnalysisToolKit
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

# Filer: B6N == B6N && BFMI == BFMI (not on X, Y)
# Filter first: B6N != BFMI
# Then: Filter 3-alleles (perhaps keep them in as 2 alleles when 1 is very low < 2 %)
# Then: Find the Frequencies of the remaining SNPs in F1 (by group)

createNames <- function(x){ paste0(x[,1],":", x[,2],"_", x[,5]) }

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo   <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
mlength   <- max(chrInfo[,"Length"])
chromosomes  <- as.character(c(1:19, "X", "Y", "MT"))

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")

B6Nm <- read.table("Analysis/5068_GAGTGG_L004_.snps.vcf", colClasses="character")
B6Nf <- read.table("Analysis/5069_AGTCAA_L004_.snps.vcf", colClasses="character")

namesB6Nf <- createNames(B6Nf) ; namesB6Nm <- createNames(B6Nm)

B6N <- B6Nm[which(namesB6Nm %in% namesB6Nf),]
colnames(B6N) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

BFMIm <- read.table("Analysis/4868_GCCAAT_L001_.snps.vcf", colClasses="character")
BFMIf <- read.table("Analysis/5067_ATCACG_L004_.snps.vcf", colClasses="character")

namesBFMIf <- createNames(BFMIf) ; namesBFMIm <- createNames(BFMIm)

BFMI <- BFMIm[which(namesBFMIm %in% namesBFMIf),]
colnames(BFMI) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

# maternal B6N

matB6_1 <- read.table("Analysis/5070_CGATGT_L005_.snps.vcf", colClasses="character")
matB6_2 <- read.table("Analysis/5071_CCGTCC_L005_.snps.vcf", colClasses="character")
matB6_3 <- read.table("Analysis/5072_TAGCTT_L005_.snps.vcf", colClasses="character")

namesmB6_1 <- createNames(matB6_1) ; namesmB6_2 <- createNames(matB6_2) ; namesmB6_3 <- createNames(matB6_3)
matB6 <- matB6_1[which(namesmB6_1 %in% namesmB6_2 & namesmB6_1 %in% namesmB6_3),]
colnames(matB6) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

# maternal BFMI

matBFMI_1 <- read.table("Analysis/5073_TTAGGC_L006_.snps.vcf", colClasses="character")
matBFMI_2 <- read.table("Analysis/5074_GATCAG_L006_.snps.vcf", colClasses="character")
matBFMI_3 <- read.table("Analysis/5075_ATGTCA_L006_.snps.vcf", colClasses="character")

namesmBFMI_1 <- createNames(matBFMI_1) ; namesmBFMI_2 <- createNames(matBFMI_2) ; namesmBFMI_3 <- createNames(matBFMI_3)
matBFMI <- matBFMI_1[which(namesmBFMI_1 %in% namesmBFMI_2 & namesmBFMI_1 %in% namesmBFMI_3),]
colnames(matBFMI) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

rownames(B6N)   <- createNames(B6N)
rownames(BFMI)   <- createNames(BFMI)
rownames(matB6)   <- createNames(matB6)
rownames(matBFMI) <- createNames(matBFMI)

doAnalysis <- function(maternal, BFMI, B6N){
  mmatrix <- NULL
  for(snp in rownames(maternal)){
    if(maternal[snp,"FORMAT"] == "GT:AD:DP:GQ:PL"){
      values <- strsplit(maternal[snp,"SAMPLE"], ":")
      totalReads <- as.numeric(unlist(values[[1]][3]))
      if(totalReads > 10){                                                                                  # Minimum of 10 reads (combined for the alleles)
        alleleReads <- as.numeric(unlist(strsplit(unlist(values[[1]][2]),",")))
        if(length(alleleReads) == 2){                                                                       # Limit to bi-allelic SNPs
          inBFMI  <- which(rownames(BFMI) ==  snp)
          inB6N   <- which(rownames(B6N) ==  snp)
          if(length(inBFMI) > 0 && length(inB6N) == 0){                                                     # SNP found in BFMI, not B6N
            cat(snp,": ", totalReads,"->", alleleReads/totalReads, "\n")
            mmatrix <- rbind(mmatrix, c(snp, maternal[snp,"CHROM"], maternal[snp,"POS"], maternal[snp,"ID"], "BFMI", inBFMI, alleleReads/totalReads))
          }
          if(length(inBFMI) == 0 && length(inB6N) > 0){                                                     # SNP found in B6N, not BFMI
            cat(snp,": ", totalReads,"->", alleleReads/totalReads, "\n")
            mmatrix <- rbind(mmatrix, c(snp, maternal[snp,"CHROM"], maternal[snp,"POS"], maternal[snp,"ID"], "B6N", inB6N, alleleReads/totalReads))
          }
        }
      }
    }
   # if(!is.null(mmatrix) && nrow(mmatrix) > 1000){
   #   colnames(mmatrix) <- c("ID", "Chr", "Loc", "dbSNP", "Origin", "OriginLoc", "Reference", "Alternative")
   #  return(mmatrix)
   # }
  }
  colnames(mmatrix) <- c("ID", "Chr", "Loc", "dbSNP", "Origin", "OriginLoc", "Reference", "Alternative")
  return(mmatrix)
}

matB6Nsnps <- doAnalysis(matB6, BFMI, B6N)
matBFMIsnps <- doAnalysis(matBFMI, BFMI, B6N)

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="SNP origin", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1

aa <- apply(matB6Nsnps, 1,function(x){
  yloc <- match(as.character(x["Chr"]), chromosomes); xloc <- as.numeric(x["Loc"])
  col <- "white"
  if(as.numeric(x["Alternative"]) > 0.9 && x["Origin"] == "BFMI") col <- "orange"
  if(as.numeric(x["Alternative"]) > 0.9 && x["Origin"] == "B6N") col <- "gray"
  if(col != "white") points(x=xloc, y=yloc - 0.1, pch=15,cex=0.5, col=col)
})

aa <- apply(matBFMIsnps, 1,function(x){
  yloc <- match(as.character(x["Chr"]), chromosomes); xloc <- as.numeric(x["Loc"])
  col <- "white"
  if(as.numeric(x["Alternative"]) > 0.9 && x["Origin"] == "BFMI") col <- "orange"
  if(as.numeric(x["Alternative"]) > 0.9 && x["Origin"] == "B6N") col <- "gray"
  if(col != "white") points(x=xloc, y=yloc + 0.1, pch=15,cex=0.5, col=col)
})

aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1,lwd=2)
  cnt <<- cnt + 1
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
legend("topright", c("> 90% BFMI", "> 90% B6N"), fill=c("orange","gray"))

ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

chr <- 6
op <- par(mfrow=c(2,1))
mBFMI <- matBFMIsnps[matBFMIsnps[,"Chr"]==chr,"Alternative"]
plot(mBFMI, col=as.numeric(as.factor(matBFMIsnps[,"Origin"])), pch=19,cex=1)
points(ma(mBFMI, 25), t='l',lwd=2)

mB6N <- matB6Nsnps[matB6Nsnps[,"Chr"]==chr,"Alternative"]
plot(mB6N, col=as.numeric(as.factor(matB6Nsnps[,"Origin"])), pch=19,cex=1)
points(ma(mB6N, 25), t='l',lwd=2)

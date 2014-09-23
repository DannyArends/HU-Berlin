# snpanalysis.R - Analyze the SNPs and indels called by the GenomeAnalysisToolKit
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Sep, 2014

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

colnames(B6Nm) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(B6Nf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

BFMIm <- read.table("Analysis/4868_GCCAAT_L001_.snps.vcf", colClasses="character")
BFMIf <- read.table("Analysis/5067_ATCACG_L004_.snps.vcf", colClasses="character")

colnames(BFMIm) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(BFMIf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

matB6_1 <- read.table("Analysis/5070_CGATGT_L005_.snps.vcf", colClasses="character")        # maternal B6N
matB6_2 <- read.table("Analysis/5071_CCGTCC_L005_.snps.vcf", colClasses="character")        # maternal B6N
matB6_3 <- read.table("Analysis/5072_TAGCTT_L005_.snps.vcf", colClasses="character")        # maternal B6N

rownames(matB6_1) <- createNames(matB6_1) ; rownames(matB6_2) <- createNames(matB6_2) ; rownames(matB6_3) <- createNames(matB6_3)
matB6 <- matB6_1[which(rownames(matB6_1) %in% rownames(matB6_2) & rownames(matB6_1) %in% rownames(matB6_3)), ]
colnames(matB6) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_1) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_3) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

matBFMI_1 <- read.table("Analysis/5073_TTAGGC_L006_.snps.vcf", colClasses="character")      # maternal BFMI
matBFMI_2 <- read.table("Analysis/5074_GATCAG_L006_.snps.vcf", colClasses="character")      # maternal BFMI
matBFMI_3 <- read.table("Analysis/5075_ATGTCA_L006_.snps.vcf", colClasses="character")      # maternal BFMI

rownames(matBFMI_1) <- createNames(matBFMI_1) ; rownames(matBFMI_2) <- createNames(matBFMI_2) ; rownames(matBFMI_3) <- createNames(matBFMI_3)
matBFMI <- matBFMI_1[which(rownames(matBFMI_1) %in% rownames(matBFMI_2) & rownames(matBFMI_1) %in% rownames(matBFMI_3)), ]
colnames(matBFMI) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_1) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_3) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

rownames(B6Nm)     <- createNames(B6Nm)   ; rownames(B6Nf)     <- createNames(B6Nf)
rownames(BFMIm)    <- createNames(BFMIm)  ; rownames(BFMIf)    <- createNames(BFMIf)
rownames(matB6)    <- createNames(matB6)  ; rownames(matBFMI)  <- createNames(matBFMI)

doAnalysis <- function(maternal, m1, m2, m3){
  mmatrix <- NULL
  for(snp in rownames(maternal)){
    if(maternal[snp,"FORMAT"] == "GT:AD:DP:GQ:PL"){
      v1 <- strsplit(m1[snp,"SAMPLE"], ":") ; v2 <- strsplit(m2[snp,"SAMPLE"], ":") ; v3 <- strsplit(m3[snp,"SAMPLE"], ":")
      v1Reads <- as.numeric(unlist(v1[[1]][3])) ; v2Reads <- as.numeric(unlist(v2[[1]][3])) ; v3Reads <- as.numeric(unlist(v3[[1]][3]))
      if(v1Reads > 5 && v2Reads > 5 && v3Reads > 5){                                                        # Minimum of 5 reads (combined for the alleles)
        v1ReadsA <- as.numeric(unlist(strsplit(unlist(v1[[1]][2]),",")))
        v2ReadsA <- as.numeric(unlist(strsplit(unlist(v2[[1]][2]),",")))
        v3ReadsA <- as.numeric(unlist(strsplit(unlist(v3[[1]][2]),",")))
        if(length(v1ReadsA) == 2 && length(v2ReadsA) == 2 && length(v3ReadsA) == 2){                        # Limit to bi-allelic SNPs
          inBFMI  <- c(which(rownames(BFMIm) ==  snp), which(rownames(BFMIf) ==  snp))                      # SNP in BFMI males/females
          inB6N   <- c(which(rownames(B6Nm) ==  snp), which(rownames(B6Nf) ==  snp))                        # SNP in B6N males/females
          if(length(inBFMI) == 2 && length(inB6N) == 0){                                                    # SNP found in BFMI, not B6N
            impScore <- (abs((v1ReadsA/v1Reads)[2] - 0.5) + abs((v2ReadsA/v2Reads)[2] - 0.5) + abs((v3ReadsA/v3Reads)[2] - 0.5)) / 3
            cat(snp,": ", v1Reads, v2Reads, v3Reads,"->", v1ReadsA/v1Reads, v2ReadsA/v2Reads, v3ReadsA/v3Reads, ":", impScore, "\n")
            mmatrix <- rbind(mmatrix, c(snp, maternal[snp,"CHROM"], maternal[snp,"POS"], maternal[snp,"ID"], "BFMI", inBFMI, v1ReadsA/v1Reads, v2ReadsA/v2Reads, v3ReadsA/v3Reads, impScore))
          }
          if(length(inBFMI) == 0 && length(inB6N) == 2){                                                    # SNP found in B6N, not BFMI
            impScore <- (abs((v1ReadsA/v1Reads)[2] - 0.5) + abs((v2ReadsA/v2Reads)[2] - 0.5) + abs((v3ReadsA/v3Reads)[2] - 0.5)) / 3
            cat(snp,": ", v1Reads, v2Reads, v3Reads,"->", v1ReadsA/v1Reads, v2ReadsA/v2Reads, v3ReadsA/v3Reads, ":", impScore, "\n")
            mmatrix <- rbind(mmatrix, c(snp, maternal[snp,"CHROM"], maternal[snp,"POS"], maternal[snp,"ID"], "B6N",  inB6N, v1ReadsA/v1Reads, v2ReadsA/v2Reads, v3ReadsA/v3Reads, impScore))
          }
        }
      }
    }
    #if(!is.null(mmatrix) && nrow(mmatrix) > 5000){
    #  colnames(mmatrix) <- c("ID", "Chr", "Loc", "dbSNP", "Origin", "OriginPaternal", "OriginMaternal", "R1", "A1", "R2", "A2", "R3", "A3", "ImprintingScore")
    #  return(mmatrix)
    #}
  }
  colnames(mmatrix) <- c("ID", "Chr", "Loc", "dbSNP", "Origin", "OriginPaternal", "OriginMaternal", "R1", "A1", "R2", "A2", "R3", "A3",  "ImprintingScore")
  return(mmatrix)
}

matB6Nsnps <- doAnalysis(matB6, matB6_1, matB6_2, matB6_3)
matBFMIsnps <- doAnalysis(matBFMI, matBFMI_1, matBFMI_2, matBFMI_3)

write.table(matB6Nsnps, file="maternalB6snps_5reads.txt", sep="\t", row.names=FALSE)                                  # Also available for 10 reads per individual
write.table(matBFMIsnps, file="maternalBFMIsnps_5reads.txt", sep="\t", row.names=FALSE)                               # Also available for 10 reads per individual

### Downstream analysis and some plots
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
matB6Nsnps  <- read.table("maternalB6snps_5reads.txt", sep="\t", header=TRUE)
matBFMIsnps <- read.table("maternalBFMIsnps_5reads.txt", sep="\t", header=TRUE)

### PLOTS

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="SNP origin", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1

aa <- apply(matB6Nsnps, 1,function(x){
  yloc <- match(as.character(x["Chr"]), chromosomes); xloc <- as.numeric(x["Loc"])
  col <- "white"
  if(as.numeric(x["ImprintingScore"]) > 0.3 && x["Origin"] == "BFMI") col <- "orange"
  if(as.numeric(x["ImprintingScore"]) > 0.3 && x["Origin"] == "B6N") col <- "gray"
  if(col != "white") points(x=xloc, y=yloc - 0.1, pch=15,cex=0.5, col=col)
})

aa <- apply(matBFMIsnps, 1,function(x){
  yloc <- match(as.character(x["Chr"]), chromosomes); xloc <- as.numeric(x["Loc"])
  col <- "white"
  if(as.numeric(x["ImprintingScore"]) > 0.3 && x["Origin"] == "BFMI") col <- "orange"
  if(as.numeric(x["ImprintingScore"]) > 0.3 && x["Origin"] == "B6N") col <- "gray"
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


plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="SNP origin", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1

aa <- apply(matBFMIsnps, 1,function(x){
  yloc <- match(as.character(x["Chr"]), chromosomes); xloc <- as.numeric(x["Loc"])
  col <- "black"
  if(as.numeric(x["ImprintingScore"]) > 0.1 && x["Origin"] == "BFMI") col <- "orange"
  if(as.numeric(x["ImprintingScore"]) > 0.1 && x["Origin"] == "B6N") col <- "gray"
  points(x=xloc, y=yloc+as.numeric(x["ImprintingScore"]), pch=19, cex=0.3, col=col)
})

aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1,lwd=2)
  cnt <<- cnt + 1
})

mB6N <- matB6Nsnps[matB6Nsnps[,"Chr"]==chr,"ImprintingScore"]
plot(mB6N, col=as.numeric(as.factor(matB6Nsnps[,"Origin"])), pch=19,cex=1)
points(ma(mB6N, 25), t='l',lwd=2)

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

B6Nm <- read.table("Analysis/5068_GAGTGG_L004_.snps.bcftools.vcf", colClasses="character")
B6Nf <- read.table("Analysis/5069_AGTCAA_L004_.snps.bcftools.vcf", colClasses="character")

colnames(B6Nm) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(B6Nf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

BFMIm <- read.table("Analysis/4868_GCCAAT_L001_.snps.bcftools.vcf", colClasses="character")
BFMIf <- read.table("Analysis/5067_ATCACG_L004_.snps.bcftools.vcf", colClasses="character")

colnames(BFMIm) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(BFMIf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

matB6_1 <- read.table("Analysis/5070_CGATGT_L005_.snps.bcftools.vcf", colClasses="character")        # maternal B6N
matB6_2 <- read.table("Analysis/5071_CCGTCC_L005_.snps.bcftools.vcf", colClasses="character")        # maternal B6N
matB6_3 <- read.table("Analysis/5072_TAGCTT_L005_.snps.bcftools.vcf", colClasses="character")        # maternal B6N

rownames(matB6_1) <- createNames(matB6_1) ; rownames(matB6_2) <- createNames(matB6_2) ; rownames(matB6_3) <- createNames(matB6_3)
matB6 <- matB6_1[which(rownames(matB6_1) %in% rownames(matB6_2) & rownames(matB6_1) %in% rownames(matB6_3)), ]
colnames(matB6) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_1) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_3) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

matBFMI_1 <- read.table("Analysis/5073_TTAGGC_L006_.snps.bcftools.vcf", colClasses="character")      # maternal BFMI
matBFMI_2 <- read.table("Analysis/5074_GATCAG_L006_.snps.bcftools.vcf", colClasses="character")      # maternal BFMI
matBFMI_3 <- read.table("Analysis/5075_ATGTCA_L006_.snps.bcftools.vcf", colClasses="character")      # maternal BFMI

rownames(matBFMI_1) <- createNames(matBFMI_1) ; rownames(matBFMI_2) <- createNames(matBFMI_2) ; rownames(matBFMI_3) <- createNames(matBFMI_3)
matBFMI <- matBFMI_1[which(rownames(matBFMI_1) %in% rownames(matBFMI_2) & rownames(matBFMI_1) %in% rownames(matBFMI_3)), ]
colnames(matBFMI) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_1) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_3) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

rownames(B6Nm)     <- createNames(B6Nm)   ; rownames(B6Nf)     <- createNames(B6Nf)
rownames(BFMIm)    <- createNames(BFMIm)  ; rownames(BFMIf)    <- createNames(BFMIf)
rownames(matB6)    <- createNames(matB6)  ; rownames(matBFMI)  <- createNames(matBFMI)

cat("SNPs matBFMI",nrow(matBFMI_1),nrow(matBFMI_2),nrow(matBFMI_3), nrow(matBFMI),"\n")
cat("SNPs matBFMI on X/Y/MT:", length(c(which(matBFMI[,"CHROM"] == "X"), which(matBFMI[,"CHROM"] == "MT"), which(matBFMI[,"CHROM"] == "Y"))),"\n")
cat("SNPs matBFMI",nrow(matB6_1),nrow(matB6_2),nrow(matB6_3), nrow(matB6),"\n")
cat("SNPs matB6N on X/Y/MT:", length(c(which(matB6[,"CHROM"] == "X"), which(matB6[,"CHROM"] == "MT"), which(matB6[,"CHROM"] == "Y"))),"\n")

doAnalysis <- function(maternal, m1, m2, m3){
  mmatrix <- NULL
  for(snp in rownames(maternal)){
    if(maternal[snp,"FORMAT"] == "GT:PL"){
      v1 <- strsplit(m1[snp,"INFO"], ";") ; v2 <- strsplit(m2[snp,"INFO"], ";") ; v3 <- strsplit(m3[snp,"INFO"], ";")
      v1Reads <- as.numeric(unlist(strsplit(gsub("DP4=","",unlist(v1)[which(grepl("DP4", unlist(v1)))]),",")))
      v2Reads <- as.numeric(unlist(strsplit(gsub("DP4=","",unlist(v2)[which(grepl("DP4", unlist(v2)))]),",")))
      v3Reads <- as.numeric(unlist(strsplit(gsub("DP4=","",unlist(v3)[which(grepl("DP4", unlist(v3)))]),",")))
      if(sum(v1Reads) >= 10 && sum(v2Reads) >= 10 && sum(v3Reads) >= 10){                                                     # Minimum of 10 reads (combined for the alleles)
        v1ReadsA <- sum(v1Reads[3:4]) ; v1Reads <- sum(v1Reads[1:4])
        v2ReadsA <- sum(v2Reads[3:4]) ; v2Reads <- sum(v2Reads[1:4])
        v3ReadsA <- sum(v3Reads[3:4]) ; v3Reads <- sum(v3Reads[1:4])
      #  if(length(v1ReadsA) == 2 && length(v2ReadsA) == 2 && length(v3ReadsA) == 2){                        # Limit to bi-allelic SNPs
          inBFMI  <- c(which(rownames(BFMIm) ==  snp), which(rownames(BFMIf) ==  snp))                      # SNP in BFMI males/females
          inB6N   <- c(which(rownames(B6Nm) ==  snp), which(rownames(B6Nf) ==  snp))                        # SNP in B6N males/females
          if(length(inBFMI) == 2 && length(inB6N) == 0){                                                    # SNP found in BFMI, not B6N
            #cat("SNP in BFMI\n")
            r1 <- v1ReadsA/v1Reads; r2 <- v2ReadsA/v2Reads; r3 <- v3ReadsA/v3Reads
            impScore <- (abs(r1 - 0.5) + abs(r2 - 0.5) + abs(r3 - 0.5)) / 3
            cat(snp,": ", v1Reads, v2Reads, v3Reads,"->", r1, r2, r3, ":", impScore, "\n")
            origin <- c("BFMI", "BFMI", "BFMI")
            if(r1 < 0.5) origin[1] <- "B6N"
            if(r2 < 0.5) origin[2] <- "B6N"
            if(r3 < 0.5) origin[3] <- "B6N"
            mmatrix <- rbind(mmatrix, c(snp, maternal[snp,"CHROM"], maternal[snp,"POS"], maternal[snp,"ID"], origin, inBFMI, v1ReadsA/v1Reads, v1Reads, v2ReadsA/v2Reads, v2Reads, v3ReadsA/v3Reads, v3Reads, impScore, "BFMIsnp"))
          }
          if(length(inBFMI) == 0 && length(inB6N) == 2){                                                    # SNP found in B6N, not BFMI
            #cat("SNP in B6N\n")
            r1 <- v1ReadsA/v1Reads; r2 <- v2ReadsA/v2Reads; r3 <- v3ReadsA/v3Reads
            impScore <- (abs(r1 - 0.5) + abs(r2 - 0.5) + abs(r3 - 0.5)) / 3
            cat(snp,": ", v1Reads, v2Reads, v3Reads,"->", r1, r2, r3, ":", impScore, "\n")
            origin <- c("B6N", "B6N", "B6N")
            if(r1 < 0.5) origin[1] <- "BFMI"
            if(r2 < 0.5) origin[2] <- "BFMI"
            if(r3 < 0.5) origin[3] <- "BFMI"
            mmatrix <- rbind(mmatrix, c(snp, maternal[snp,"CHROM"], maternal[snp,"POS"], maternal[snp,"ID"], origin,  inB6N, v1ReadsA/v1Reads, v1Reads, v2ReadsA/v2Reads, v2Reads, v3ReadsA/v3Reads, v3Reads, impScore, "B6Nsnp"))
          }
      #  }
      }else{
        cat("reads FAILED",v1Reads,v2Reads,v3Reads,"\n")
      }
    }
  }
  colnames(mmatrix) <- c("ID", "Chr", "Loc", "dbSNP", "Origin1", "Origin2", "Origin3", "OriginPaternal", "OriginMaternal", "R1", "N1", "R2", "N2", "R3", "N3", "ImprintingScore", "Detected")
  return(mmatrix)
}

matB6Nsnps <- doAnalysis(matB6, matB6_1, matB6_2, matB6_3)
matBFMIsnps <- doAnalysis(matBFMI, matBFMI_1, matBFMI_2, matBFMI_3)

write.table(matB6Nsnps, file="maternalB6snps_10reads.txt", sep="\t", row.names=FALSE)                                  # Also available for 10 reads per individual
write.table(matBFMIsnps, file="maternalBFMIsnps_10reads.txt", sep="\t", row.names=FALSE)                               # Also available for 10 reads per individual

### PLOT
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
matB6Nsnps  <- read.table("maternalB6snps_5reads.txt", sep="\t", header=TRUE)
matBFMIsnps <- read.table("maternalBFMIsnps_5reads.txt", sep="\t", header=TRUE)

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

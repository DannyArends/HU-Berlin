# SNP calling using the population VCF (only heterozygous)
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

getDP <- function(x){
  v1 <- strsplit(x, ";")
  as.numeric(unlist(strsplit(gsub("DP=","",unlist(v1)[grep("DP=", unlist(v1))]),",")))
}

createNames <- function(x){ paste0(x[,1],":", x[,2],"_", x[,5]) }
getGenotypes <- function(x){ unlist(lapply(strsplit(x, ":"),"[",1)) }
getProbabilities <- function(x){ unlist(lapply(strsplit(x, ":"),"[",2)) }
getMaxProb <- function(x){ 
  p <- getProbabilities(x)
  unlist(lapply(strsplit(p, ","),function(i){max(as.numeric(i))}))
}

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
vcfdata <- read.table("population.vcf", colClasses="character")
colnames(vcfdata) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","5068","5069","4868","5067","5070","5071","5072","5073","5074","5075")
rownames(vcfdata) <- createNames(vcfdata)
parents <- c("5068","5069","4868","5067")
matB6N  <- c("5070","5071","5072")
matBFMI <- c("5073","5074","5075")
samples <- c("5068","5069","4868","5067","5070","5071","5072","5073","5074","5075")

DPs <- unlist(lapply(vcfdata[,"INFO"], getDP))
vcfdata <- vcfdata[DPs > 100,]

genos <- matrix(NA, nrow(vcfdata), length(samples))
probs <- matrix(NA, nrow(vcfdata), length(samples))
colnames(genos) <- samples ; rownames(genos) <- createNames(vcfdata)
colnames(probs) <- samples ; rownames(probs) <- createNames(vcfdata)

for(s in samples){ 
  genos[,s] <- getGenotypes(vcfdata[,s])
  probs[,s] <- getMaxProb(vcfdata[,s])
}

confidenceThreshold <- 50

probInParents <- which(apply(probs[,parents], 1,function(x){ sum(x > confidenceThreshold) == 4 }))
genos <- genos[probInParents, ] ; probs <- probs[probInParents, ]

difGenoInParents <- which(apply(genos[,parents], 1,function(x){ 
  if(any(x == "0/1")) return(FALSE)
  if(x[1] == x[2] && x[3] == x[4] && x[1] != x[3]) return(TRUE)
  return(FALSE)
}))

genos <- genos[difGenoInParents, ] ; probs <- probs[difGenoInParents, ]
genos[which(probs < confidenceThreshold)] <- NA

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
matB6_1 <- read.table("Analysis/5070_CGATGT_L005_.snps.bcftools.vcf", colClasses="character")        # maternal B6N
matB6_2 <- read.table("Analysis/5071_CCGTCC_L005_.snps.bcftools.vcf", colClasses="character")        # maternal B6N
matB6_3 <- read.table("Analysis/5072_TAGCTT_L005_.snps.bcftools.vcf", colClasses="character")        # maternal B6N
rownames(matB6_1) <- createNames(matB6_1) ; rownames(matB6_2) <- createNames(matB6_2) ; rownames(matB6_3) <- createNames(matB6_3)
matB6 <- matB6_1[which(rownames(matB6_1) %in% rownames(matB6_2) & rownames(matB6_1) %in% rownames(matB6_3)), ]

colnames(matB6_1) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matB6_3) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

matBFMI_1 <- read.table("Analysis/5073_TTAGGC_L006_.snps.bcftools.vcf", colClasses="character")      # maternal BFMI
matBFMI_2 <- read.table("Analysis/5074_GATCAG_L006_.snps.bcftools.vcf", colClasses="character")      # maternal BFMI
matBFMI_3 <- read.table("Analysis/5075_ATGTCA_L006_.snps.bcftools.vcf", colClasses="character")      # maternal BFMI
rownames(matBFMI_1) <- createNames(matBFMI_1) ; rownames(matBFMI_2) <- createNames(matBFMI_2) ; rownames(matBFMI_3) <- createNames(matBFMI_3)
matBFMI <- matBFMI_1[which(rownames(matBFMI_1) %in% rownames(matBFMI_2) & rownames(matBFMI_1) %in% rownames(matBFMI_3)), ]

colnames(matBFMI_1) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
colnames(matBFMI_3) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

getRatio <- function(X, genos){
  X_1 <- strsplit(X[rownames(genos)[x],"INFO"], ";")
  X_1r <- as.numeric(unlist(strsplit(gsub("DP4=","",unlist(X_1)[which(grepl("DP4", unlist(X_1)))]),",")))
  X_1alt <- sum(X_1r[3:4]) ; X_1all <- sum(X_1r[1:4]) ; X_1rat <- X_1alt / X_1all
  return(X_1rat)
}

SNPs <- NULL
cnt <- 0
for(x in 1:nrow(genos)){
  inB6N <- rownames(genos)[x] %in% rownames(matB6)
  inBFMI <- rownames(genos)[x] %in% rownames(matBFMI)
  if(inB6N && inBFMI){
    b6rat <- c(getRatio(matB6_1,genos), getRatio(matB6_2,genos), getRatio(matB6_3,genos))
    bfrat <- c(getRatio(matBFMI_1,genos), getRatio(matBFMI_2,genos), getRatio(matBFMI_3,genos))
    if(sum(b6rat < 0.5) == 3 || sum(b6rat > 0.5) == 3){
      if(sum(bfrat < 0.5) == 3 || sum(bfrat > 0.5) == 3){
        if(abs(mean(b6rat) - mean(bfrat)) > 0.2){
          cat(rownames(genos)[x], mean(b6rat), mean(bfrat), "\n")
          cnt <- cnt + 1
          SNPs <- rbind(SNPs, c(rownames(genos)[x], matB6_1[rownames(genos)[x],"CHROM"], matB6_1[rownames(genos)[x],"POS"], mean(b6rat), mean(bfrat)))
        }
      }
    }
  }
}
colnames(SNPs) <- c("snpID", "CHROM", "POS", "matB6Nratio", "matBFMIratio")

geneSNPCoupling <- vector("list", length(uniqueGenes))
x <- 1
for(gene in geneExonsCoupling){                                                                               # Summarize all SNPs from a direction per gene
  onChr <- as.character(SNPs[,"CHROM"]) == as.character(unique(gene[,3]))
  for(exon in 1:nrow(gene)){
    exonStart <- as.numeric(gene[exon,6])
    exonEnd <- as.numeric(gene[exon,7])
    inEXON <- which(SNPs[onChr,"POS"] >= exonStart & SNPs[onChr,"POS"] <= exonEnd)
    if(length(inEXON) > 0){
      geneSNPCoupling[[x]] <- rbind(geneSNPCoupling[[x]], as.matrix(SNPs[onChr,])[inEXON,])
      geneSNPCoupling[[x]] <- geneSNPCoupling[[x]][!duplicated(geneSNPCoupling[[x]][,"snpID"]),]
    }
    #names(geneSNPCoupling)[x] <- as.character(unique(gene[,1]))
  }
  
  if(!is.null(geneSNPCoupling[[x]])){
    if(is.null(nrow(geneSNPCoupling[[x]]))) geneSNPCoupling[[x]] <- t(as.matrix(geneSNPCoupling[[x]]))
    cat(x, as.character(unique(gene[,1])), "found", nrow(geneSNPCoupling[[x]]), ncol(geneSNPCoupling[[x]]), "SNPs in gene\n")
  }
  x <- x + 1
}

names(geneSNPCoupling) <- uniqueGenes

geneSNPCoupling$ENSMUSG00000040102 <- t(as.matrix(SNPs[18,]))
geneSNPCoupling$ENSMUSG00000024370 <- t(as.matrix(SNPs[54,]))

geneSNPshort <- NULL
for(x in 1:length(geneSNPCoupling)){
  if(!is.null(geneSNPCoupling[[x]])) geneSNPshort <- c(geneSNPshort, geneSNPCoupling[x])
}

matrixform <- NULL
for(x in 1:length(geneSNPshort)){
  additionalgeneInfo <- RPKM[which(RPKM[,"ensembl_gene_id"] == names(geneSNPshort)[x]), -c(1)]
  matrixform <- rbind(matrixform, cbind(ensembl_gene_id = names(geneSNPshort)[x], cbind(geneSNPshort[[x]]), additionalgeneInfo))
}

write.table(matrixform,"RPKM+ASE_RATIOS_MinReads100_Threshold50.txt", sep="\t", quote=FALSE,row.names=FALSE)
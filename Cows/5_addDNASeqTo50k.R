
# Write out the file to scale down the complete SNP calling
write.table(cbind(as.character(snpinfo[,"Chr"]), snpinfo[,"SNPpos"]), "vcf_filter.txt",sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
# Copy to server, run vcftools to extract a subset
#vcftools --vcf /home/danny/NAS/Cattle/DSN/SNPcalling/combined/population41/0000_gatk.vcf --positions vcf_filter.txt --recode --out subset_0000_gatk.vcf
# Copy results back to my pc, remove the # in front of the header

## Merge with DNA-Seq data
setwd("D:/Ddrive/GiesenFugato")
all50Kdata <- read.table("merged_GVF.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
snpinfo <- read.table("merged_SNPinfo.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)

seqVCF <- read.table("DNASeq/subset_0000_gatk.vcf.recode.vcf", header=TRUE, colClasses="character", stringsAsFactor=FALSE)
complexAlleles <-  grep(",", seqVCF[,"ALT"])
seqVCF <- seqVCF[-complexAlleles, ]
seqVCF <- cbind(ChrPos = paste(seqVCF[,"CHROM"], seqVCF[,"POS"], sep="_"), seqVCF)
write.table(seqVCF, "DNASeq/subset_0000_gatk.vcf.recode.matched.txt", sep="\t", quote=FALSE)

matchedSNPs <- NULL
for(x in 1:nrow(snpinfo)){
  ChrPos <- paste(snpinfo[x,"Chr"], snpinfo[x,"SNPpos"], sep="_")
  inVCF <- which(seqVCF[, "ChrPos"] == ChrPos)
  if(length(inVCF) == 1) {
    matchedSNPs <- c(matchedSNPs, rownames(snpinfo)[x])
  }
}

matched50Kdata <- all50Kdata[matchedSNPs,]
matchedSNPinfo <- snpinfo[matchedSNPs,]

write.table(matched50Kdata, file = "matched_GVF_DNAseq.txt", sep="\t", quote=FALSE)
write.table(matchedSNPinfo, file = "matched_SNPinfo_DNAseq.txt", sep="\t", quote=FALSE)

## Combine / Merge all data in VCF coding 
opposite <- function(x){
  if(x == "T") return("A")
  if(x == "C") return("G")
  if(x == "A") return("T")
  if(x == "G") return("C")
}

setwd("D:/Ddrive/GiesenFugato")
matched50Kdata <- read.table("matched_GVF_DNAseq.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
matchedSNPinfo <- read.table("matched_SNPinfo_DNAseq.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
seqVCF <- read.table("DNASeq/subset_0000_gatk.vcf.recode.matched.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)

nsamples <- ncol(matched50Kdata) + (ncol(seqVCF)-10)

mergedData <- matrix(NA, nrow(matched50Kdata), nsamples, dimnames=list(rownames(matched50Kdata), c(colnames(matched50Kdata), colnames(seqVCF)[-c(1:10)])))

for(x in 1:nrow(matchedSNPinfo)){
  ChrPos <- paste(matchedSNPinfo[x,"Chr"], matchedSNPinfo[x,"SNPpos"], sep="_")
  inVCF <- which(seqVCF[, "ChrPos"] == ChrPos)
  if(length(inVCF) == 1) {
    GTsequencing <- unlist(lapply(strsplit(as.character(seqVCF[inVCF, -c(1:10)]), ":"),"[",1))
    ref <- seqVCF[inVCF,"REF"]
    alt <- seqVCF[inVCF,"ALT"]
    
    GT50K <- strsplit(as.character(matched50Kdata[x, ]), "")
    snpcodes <- names(table(unlist(GT50K)))
    
    if(!(ref %in% snpcodes)) {    # Reference allele not in the encoding of the SNP genotype, Flip the reference allele
      ref <- opposite(ref)
    }
    if(!(alt %in% snpcodes)) {    # Alternative allele not in the encoding of the SNP genotype, Flip the alternative allele
      alt <- opposite(alt)
    }
    if(!(ref %in% snpcodes) || !(alt %in% snpcodes) || ref == alt){
      if(ref == alt){
        cat("SKIPPING ", rownames(matchedSNPinfo)[x], ", ref == alt\n")
      }else{
        cat("SNP coding ",paste0(snpcodes), "mismatches at ", rownames(matchedSNPinfo)[x], ", ref == ",ref,"\n")
        cat("SNP coding ",paste0(snpcodes), "mismatches at ", rownames(matchedSNPinfo)[x], ", alt == ",alt,"\n")
      }
    }else{
      SNpGTvcf <- unlist(lapply(GT50K, function(gts){
        if(length(gts) != 2){ return(NA); }
        if(gts[1] == ref && gts[2] == ref) return("0/0");
        if(gts[1] == ref && gts[2] == alt) return("0/1");
        if(gts[1] == alt && gts[2] == ref) return("0/1");
        if(gts[1] == alt && gts[2] == alt) return("1/1");
        cat("SHOULD NEVER BE HERE")
        return(NA);
      }))
      mergedData[x,] <- c(SNpGTvcf, GTsequencing)
    }
  }else{
    stop("Should not happen on matched data")
  }
}

write.table(mergedData, file = "merged_GVF_DNAseq.txt", sep="\t", quote=FALSE)


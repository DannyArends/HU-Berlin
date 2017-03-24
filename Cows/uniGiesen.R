#
# Figure out the Uni Giessen dataformat
#


### Load in the files from Uni Giessen
setwd("D:/Ddrive/GiesenFugato")
files <- c(paste0("776/", dir("776")), paste0("842/", dir("842")))
results <- vector("list", length(files))
for(i in 1:length(files)) {
  f <- files[i]
  dataf <- read.table(f, skip=10, sep="\t", header=TRUE)
  individualID <- as.character(unique(dataf[,2]))
  GTs <- apply(dataf[,c("Allele1...Top","Allele2...Top")], 1, function(x){
    paste0(unlist(as.character(x)),collapse="")
  })
  names(GTs) <- dataf[,"SNP.Name"]
  GTs[which(dataf[,"GT.Score"] <= 0.7)] <- "--"
  results[[i]] <- list(f, individualID, GTs)
}

indnames <- unlist(lapply(results,"[",2))
probenames <- names(results[[1]][[3]])

resultsM <- matrix(NA, length(results[[1]][[3]]), length(indnames), dimnames=list(probenames, indnames))
x <- lapply(results, function(x){
  resultsM[,x[[2]]] <<- x[[3]][probenames]
})

write.table(resultsM, file="GT_UniGiessen.txt", sep="\t", quote=FALSE)

resultsM <- read.table("GT_UniGiessen.txt", header=TRUE)

ann776 <- read.csv("Giessen/Samples_Table_LOMs_776.csv", header=TRUE, sep = ";",check.names=FALSE, colClasses='character')
ann842 <- read.csv("Giessen/Samples_Table_LOMs_842.csv", header=TRUE, sep = ";",check.names=FALSE, colClasses='character')

ann_all <- rbind(ann776, ann842[,-c(10,11)])
annotation <- ann_all[, "LOM"]
names(annotation) <- apply(ann_all[, c("Plattencode", "Positionscode")],1,paste0, collapse="_")

colnames(resultsM) <- annotation[colnames(resultsM)]
resultsM[1:10,1:10]
write.table(resultsM, file="GT_UniGiessen_LOM.txt", sep="\t", quote=FALSE)


### Load in the files from the VIT

genotypes <- read.table("VIT/gZWSgenotypenDSNbullen_1612", colClasses="character",stringsAsFactor=FALSE)
gts <- apply(genotypes, 1, function(x){
  unlist(strsplit(x[2], ""))
})
colnames(gts) <- genotypes[,1]
vit_info <- read.table("VIT/SNPNamenV2mit_extra_Info_repaired.txt", sep=";", colClasses="character",stringsAsFactor=FALSE, header=TRUE)
vit_info <- vit_info[order(as.numeric(vit_info[,"SNP.Index"])), ]
rownames(gts) <- vit_info[,"SNP.Name"]
rownames(vit_info) <- vit_info[,"SNP.Name"]

for(x in 1:nrow(gts)){
  if(vit_info[rownames(gts)[x],"ALLELEA_TOP"] != ""){
    cod2 <- paste0(vit_info[rownames(gts)[x],"ALLELEA_TOP"], vit_info[rownames(gts)[x],"ALLELEA_TOP"])
    cod1 <- paste0(sort(c(vit_info[rownames(gts)[x],"ALLELEA_TOP"], vit_info[rownames(gts)[x],"ALLELEB_TOP"])), collapse="")
    cod0 <- paste0(vit_info[rownames(gts)[x],"ALLELEB_TOP"], vit_info[rownames(gts)[x],"ALLELEB_TOP"])
    gts[x, gts[x,] == "0"] <- rep(cod0, length(which(gts[x,] == "0")))
    gts[x, gts[x,] == "1"] <- rep(cod1, length(which(gts[x,] == "1")))
    gts[x, gts[x,] == "2"] <- rep(cod2, length(which(gts[x,] == "2")))
    gts[x, gts[x,] == "9"] <- rep(NA, length(which(gts[x,] == "9")))
  }else{
    gts[x,] <- rep(NA, ncol(gts))
  }
}
write.table(gts, file="GT_VIT.txt", sep="\t", quote=FALSE)


## Fix the wrong annotation for fugato (120 animals are missing in action)
sampleInfo <- read.csv("Fugato/sampleinfoKielAna.txt",sep="\t", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
fugato <- read.table("Fugato/genotypes.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE)
write.table(sampleInfo[colnames(fugato),], file="Fugato/sampleInfo_matched.txt", sep="\t", quote=FALSE)

### MERGE: Fugato, Giessen & VIT
setwd("D:/Ddrive/GiesenFugato")
giessen <- read.table("GT_UniGiessen_LOM.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE,na.strings=c("--"))
fugato <- read.table("Fugato/genotypes.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE)
vit <- read.table("GT_VIT.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE)

snp_info <- read.table("Fugato/FilteredLocs.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE, header=TRUE)
snpall <- cbind(snp_info, SNPpos = abs(as.numeric(snp_info[,"Start"]) + as.numeric(snp_info[,"Stop"])) / 2)
snpall <- snpall[with(snpall, order(SNPpos)), ]
chrs <- paste0("Chr", c(seq(1,29), "X"))
snp_info_sorted <- NULL
for(chr in chrs){
  snp_info_sorted <- rbind(snp_info_sorted, snpall[which(snpall[,"Chr"] == chr),])
}

probeNameFunny <- grep("^BFGL-NGS", snp_info_sorted[,"decoder_ID"])
snp_info_sorted[probeNameFunny,"decoder_ID"] <- paste0("ARS-", snp_info_sorted[probeNameFunny,"decoder_ID"])

giessenNA <- which(apply(giessen,1,function(x){sum(is.na(x))}) == ncol(giessen))
if(length(giessenNA) > 0) giessen <- giessen[-giessenNA, ]

vitNA <- which(apply(vit,1,function(x){sum(is.na(x))}) == ncol(vit))
if(length(vitNA) > 0) vit <- vit[-vitNA, ]

fugatoNA <- which(apply(fugato,1,function(x){sum(is.na(x))}) == ncol(fugato))
if(length(fugatoNA) > 0) fugato <- fugato[-fugatoNA, ]

dim(giessen)
dim(vit)
dim(fugato)

rownames(giessen)[grep("^BFGL-NGS", rownames(giessen))]
rownames(vit)[grep("^BFGL-NGS", rownames(vit))]
rownames(fugato)[grep("^BFGL-NGS", rownames(fugato))]

nShared <- table(c(rownames(giessen), rownames(vit), rownames(fugato), snp_info_sorted[,"decoder_ID"]) )
inAll4 <- names(nShared[which(nShared == 4)])

giessen <- giessen[which(rownames(giessen) %in% inAll4),]
vit <- vit[which(rownames(vit) %in% inAll4),]
fugato <- fugato[which(rownames(fugato) %in% inAll4),]
snpinfo <- snp_info_sorted[which(snp_info_sorted[,"decoder_ID"] %in% inAll4),]
rownames(snpinfo) <- snpinfo[,"decoder_ID"]

dim(giessen)
dim(vit)
dim(fugato)

# Order them so that all datasets match the snpinfo
giessen <- giessen[snpinfo[,"decoder_ID"], ]
vit <- vit[snpinfo[,"decoder_ID"], ]
fugato <- fugato[snpinfo[,"decoder_ID"], ]

giessen[1:5,1:5]
vit[1:5,1:5]
fugato[1:5,1:5]
snpinfo[1:5,]

write.table(giessen, file="merged_Giessen.txt", sep="\t", quote=FALSE)
write.table(vit, file="merged_VIT.txt", sep="\t", quote=FALSE)
write.table(fugato, file="merged_Fugato.txt", sep="\t", quote=FALSE)
write.table(snpinfo[,-1], file="merged_SNPinfo.txt", sep="\t", quote=FALSE)

all50Kdata <- cbind(giessen, vit, fugato)

write.table(all50Kdata, file="merged_GVF.txt", sep="\t", quote=FALSE)


# Write out the file to scale down the complete SNP calling
write.table(cbind(as.character(snpinfo[,"Chr"]), snpinfo[,"SNPpos"]), "vcf_filter.txt",sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
# Copy to server, run vcftools to extract a subset
#vcftools --vcf /home/danny/NAS/Cattle/DSN/SNPcalling/combined/population41/population.combined.vcf --positions vcf_filter.txt --recode --out subset_population.combined.vcf
# Copy results back to my pc, remove the # in front of the header

## Merge with DNA-Seq data
setwd("D:/Ddrive/GiesenFugato")
all50Kdata <- read.table("merged_GVF.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
snpinfo <- read.table("merged_SNPinfo.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
#snpinfoPaula <- read.table("paula/snpinfo50k_UG_VIT_HU.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)

all(rownames(snpinfo) %in% snpinfoPaula[, "Name"])
snpinfoPaula[!(snpinfoPaula[, "Name"] %in% rownames(snpinfo)),]


seqVCF <- read.table("DNASeq/subset_population.combined.vcf.recode.vcf", header=TRUE, colClasses="character", stringsAsFactor=FALSE)
complexAlleles <-  grep(",", seqVCF[,"ALT"])
if(length(complexAlleles) > 0 ) seqVCF <- seqVCF[-complexAlleles, ]
seqVCF <- cbind(ChrPos = paste(seqVCF[,"CHROM"], seqVCF[,"POS"], sep="_"), seqVCF)
write.table(seqVCF, "DNASeq/subset_population.combined.vcf.recode.matched.txt", sep="\t", quote=FALSE)

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

# Do some QC checks since not all SNPs merge correctly
numNAs <- apply(mergedData,1,function(x){ return(sum(is.na(x))); })
tooMuchMissing <- which(numNAs / ncol(mergedData) >= .95)

mergedData <- mergedData[-tooMuchMissing, ]
mergedSNPinfo <- matchedSNPinfo[-tooMuchMissing,]

mergedSNPinfoExtended <- matrix(NA, nrow(mergedSNPinfo), ncol(mergedSNPinfo) + 6, dimnames=list(rownames(mergedSNPinfo), c(colnames(mergedSNPinfo), colnames(seqVCF)[1:6] )))
for(x in 1:nrow(mergedSNPinfo)){
  ChrPos <- paste(mergedSNPinfo[x,"Chr"], mergedSNPinfo[x,"SNPpos"], sep="_")
  inVCF <- which(seqVCF[, "ChrPos"] == ChrPos)
  mergedSNPinfoExtended[x,] <- unlist(c(mergedSNPinfo[x,], seqVCF[inVCF, 1:6]))
}

write.table(mergedData, file = "merged_GVF_DNAseq_afterQC.txt", sep="\t", quote=FALSE)
write.table(mergedSNPinfoExtended, file = "merged_SNPinfo_afterQC.txt", sep="\t", quote=FALSE)


##

setwd("D:/Ddrive/GiesenFugato")
mergedData <- read.table("merged_GVF_DNAseq_afterQC.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
mergedSNPinfoExtended <- read.table("merged_SNPinfo_afterQC.txt", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)

# Sample information
giessen <- read.table("mergedData_LOM.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE,na.strings=c("--"))
vit <- read.table("recoded_VIT.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE)
sampleInfo <- read.csv("Fugato/sampleInfo_matched.txt",sep="\t", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)

HF25Fugato <- sample(rownames(sampleInfo)[which(sampleInfo[,"breed"] == "HF")], 25)
DSNgiesen <- colnames(giessen)
DSNvit <- colnames(vit)
DSNfugato <- rownames(sampleInfo)[which(sampleInfo[,"breed"] == "DSN")]
DSNsequencing <- colnames(mergedData)[3888:ncol(mergedData)]


## Check in the end by Dendrogram
numData <- apply(mergedData, 1, function(x){
  as.numeric(as.factor(x))
})
rownames(numData) <- colnames(mergedData)

clGiessen <- hclust(dist(numData[c(HF25Fugato, DSNgiesen), ]))
clVIT <- hclust(dist(numData[c(HF25Fugato, DSNvit), ]))
clFUG <- hclust(dist(numData[c(HF25Fugato, DSNfugato), ]))
clSEQ <- hclust(dist(numData[c(HF25Fugato, DSNsequencing), ]))

clVITFUG <- hclust(dist(numData[c(DSNfugato, DSNvit), ]))

png("giessen.png", w = 2000, h = 800)
  plot(clGiessen, hang = -1,cex=0.7, main="DSN Giessen and 25 HF (Fugato)")
dev.off()

png("vit.png", w = 2000, h = 800)
  plot(clVIT, hang = -1,cex=0.7, main="DSN VIT and 25 HF (Fugato)")
dev.off()

png("fugato.png", w = 2000, h = 800)
  plot(clFUG, hang = -1,cex=0.7, main="DSN Fugato and 25 HF (Fugato)")
dev.off()

png("vit_fug.png", w = 2000, h = 800)
  plot(clVITFUG, hang = -1,cex=0.7, main="DSN Fugato and DSN VIT)")
dev.off()

png("sequencing.png", w = 2000, h = 800)
  plot(clSEQ, hang = -1,cex=0.7, main="DSN Sequencing and 25 HF (Fugato)")
dev.off()

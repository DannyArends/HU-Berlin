
# Do some QC checks since not all SNPs merge correctly
numNAs <- apply(mergedData,1,function(x){ return(length(which(is.na(x)))); })
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

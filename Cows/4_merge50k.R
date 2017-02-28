### MERGE: Fugato, Giessen & VIT
setwd("D:/Ddrive/GiesenFugato")
giessen <- read.table("GT_UniGiessen_LOM.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE,na.strings=c("--"))
fugato <- read.table("Fugato/genotypes.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE)
vit <- read.table("recoded_VIT.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE)

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
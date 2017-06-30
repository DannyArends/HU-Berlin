
# Merge VCF files
setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")

samples <- read.table("../sample_names.txt",header=TRUE)
samples <- read.table("../samples_description.txt",sep="\t")


CSN1S1 <- read.csv("CSN1S1.SNPs.vcf", sep="\t",comment.char="#", header=FALSE, na.string=c("./.", "."), colClasses=c("character"))
colnames(CSN1S1) <- strsplit(gsub("#","",readLines("CSN1S1.SNPs.vcf",n = 34)[34]), "\t")[[1]]

CSN2 <- read.csv("CSN2.SNPs.vcf", sep="\t",comment.char="#", header=FALSE, na.string=c("./.", "."), colClasses=c("character"))
colnames(CSN2) <- strsplit(gsub("#","",readLines("CSN2.SNPs.vcf",n = 34)[34]), "\t")[[1]]

CSN1S2 <- read.csv("CSN1S2.SNPs.vcf", sep="\t",comment.char="#", header=FALSE, na.string=c("./.", "."), colClasses=c("character"))
colnames(CSN1S2) <- strsplit(gsub("#","",readLines("CSN1S2.SNPs.vcf",n = 34)[34]), "\t")[[1]]

CSN3 <- read.csv("CSN3.SNPs.vcf", sep="\t",comment.char="#", header=FALSE, na.string=c("./.", "."), colClasses=c("character"))
colnames(CSN3) <- strsplit(gsub("#","",readLines("CSN3.SNPs.vcf",n = 34)[34]), "\t")[[1]]

STC1 <- read.csv("STC1.SNPs.vcf", sep="\t",comment.char="#", header=FALSE, na.string=c("./.", "."), colClasses=c("character"))
colnames(STC1) <- strsplit(gsub("#","",readLines("STC1.SNPs.vcf",n = 34)[34]), "\t")[[1]]

DGAT1 <- read.csv("DGAT1.SNPs.vcf", sep="\t",comment.char="#", header=FALSE, na.string=c("./.", "."), colClasses=c("character"))
colnames(DGAT1) <- strsplit(gsub("#","",readLines("DGAT1.SNPs.vcf",n = 34)[34]), "\t")[[1]]

variations <- rbind(CSN1S1, CSN2, CSN1S2, CSN3, STC1, DGAT1)
filtered <- variations[which(as.numeric(variations[,"QUAL"]) > 100),]
filtered <- filtered[-grep("INDEL", filtered[,"INFO"]),]
filtered <- filtered[-grep(",", filtered[,"ALT"]),]

cat(readLines("CSN1S1.SNPs.vcf",n=34), sep="\n", file="filteredSNPs.vcf")
write.table(filtered, file="filteredSNPs.vcf", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE, append=TRUE)

filtered[,10:ncol(filtered)] <- apply(filtered[,10:ncol(filtered)], 2, function(x){
  unlist(lapply(strsplit(x, ":"), "[", 1))
})

for(x in 1:nrow(filtered)){
  for(s in 10:ncol(filtered)){
    if(filtered[x,s] == "0/0") filtered[x,s] <- paste0(filtered[x,"REF"], filtered[x,"REF"], collapse="")
    if(filtered[x,s] == "0/1") filtered[x,s] <- paste0(sort(c(filtered[x,"REF"], filtered[x,"ALT"])), collapse="")
    if(filtered[x,s] == "1/1") filtered[x,s] <- paste0(filtered[x,"ALT"], filtered[x,"ALT"], collapse="")
    if(filtered[x,s] == "./.") filtered[x,s] <- NA
  }
}

dataperind <- apply(filtered[, 10:ncol(filtered)], 2, function(x){ sum(is.na(x)) / length(x) * 100 })
sd(diff(as.numeric(filtered[,"POS"]))[which(diff(as.numeric(filtered[,"POS"])) > 0 & diff(as.numeric(filtered[,"POS"])) < 1000)])

write.table(filtered, file="filteredSNPs.txt", sep="\t", row.names=FALSE, quote=FALSE)


### Update the GFF3 files
setwd("D:/Edrive/Goat/DNA/Sequencing")
gff3 <- read.csv("ref_ASM170441v1_top_level.gff3", sep = "\t", comment.char = "#", header = FALSE, colClasses = "character")

header <- readLines("ref_ASM170441v1_top_level.gff3",n=9)

chrN <- rbind(chr6  = c("ENA|CM004567|CM004567.1", "NC_030813.1"),
              chr8  = c("ENA|CM004569|CM004569.1", "NC_030815.1"),
              chr14 = c("ENA|CM004575|CM004575.1", "NC_030821.1"))

gff3 <- gff3[which(gff3[,1] %in% chrN[,2]),]
gff3 <- gff3[which(gff3[,3] != "region"),]
gff3 <- gff3[which(gff3[,3] != "match"),]

for(x in 1:nrow(chrN)){
  gff3[which(gff3[,1] == chrN[x,2]),1] <- chrN[x,1]
}

# Write the updated GFF3
cat(paste0(header, collapse="\n"), "\n", file = "ugff3.gff3")
write.table(gff3, file="ugff3.gff3", sep="\t", append=TRUE,col.names=FALSE, row.names=FALSE, quote=FALSE)

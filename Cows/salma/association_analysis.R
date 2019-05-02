setwd("D:/Edrive/Cow/Salma Wachow")

#Phenotype data (Obtained from Uwe)
phenotypes <- read.table("phenotypes.txt", sep="\t", header=TRUE)
phenotypes <- cbind(phenotypes, year = unlist(lapply(strsplit(as.character(phenotypes[,"birth.date"]), "/"), "[",3)))

#SNPids
SNPannotation <- read.table("SNPs.txt", sep = "\t", row.names=1)
colnames(SNPannotation) <- c("RSid", "Chr", "Pos", "A1", "A2")

#Genotypes
genotypes <- read.table("genotypes.txt", sep = "\t", header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", " "))
rownames(genotypes) <- gsub("DE", "DE00", rownames(genotypes))

rawgeno <- genotypes[, rownames(SNPannotation)]

npheno <- c("AFC", "firstMkg", "Mkg", "Mkg100", "fatkg100", "fatkg", "proteinkg1", "proteinkg", "countMast")

lactations <- unique(phenotypes[,"lactation"])

lactation <- 1
lactationdata <- phenotypes[which(phenotypes[, "lactation"] == lactation),]
rownames(lactationdata) <- lactationdata[,"id"]

nodata <- which(apply(rawgeno,2,function(x){ sum(is.na(x)) }) == nrow(rawgeno))

if(length(nodata) > 0) rawgeno <- rawgeno[,-nodata]

pvalues <- c()
for(pheno in npheno){
  pvalsPheno <- c()
  for(m in 1:ncol(rawgeno)){
    res <- anova(lm(lactationdata[rownames(rawgeno), pheno] ~ rawgeno[,m]))[[5]][1]
    pvalsPheno <- c(pvalsPheno, res)
  }
  pvalues <- rbind(pvalues, pvalsPheno)
}
rownames(pvalues) <- npheno
colnames(pvalues) <- SNPannotation[colnames(rawgeno),"RSid"]

pvalues

write.table(pvalues, file = paste0("pvalues_", lactation, ".txt"), sep = "\t", quote = FALSE)

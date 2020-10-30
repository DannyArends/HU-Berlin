library(ctl)

setwd("D:/Edrive/Mouse/S1_S2")

pheno <- read.table("allPhenotypes_final.txt", sep = "\t", row.names=1)
geno <- read.table("genotypes.cleaned.txt", sep = "\t", colClasses="character")
map <- read.table("snp_map.karl.txt", sep = ",", colClasses="character", header= TRUE, row.names=1)
map <- map[which(rownames(map) %in% rownames(geno)),]

geno <- geno[rownames(map),]

rownames(pheno) <- paste0("AIL", rownames(pheno))

pheno <- pheno[colnames(geno),]
geno <- t(geno)
genot <- apply(geno, 2, function(x){as.numeric(factor(x, levels=c("A", "H", "B")))})
rownames(genot) <- rownames(geno)


pheno <- pheno[,-which(colnames(pheno) %in% c("Sex", "ID", "WG", "Mutter", "Grandma"))]
pheno <- pheno[, c("D140","Gon","Leber")]

res <- CTLscan(genotypes=genot, phenotypes = apply(pheno,2,as.numeric), adjust = FALSE)
map[names(which(res$Leber$ctl[, "Gon"] > 4)),]

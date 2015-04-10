# Preprocessing of the MegaMuga data, creating PCA loading factors across the genotype matrix to see if they capture the population structure
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

library("psych")

setwd("E:/Mouse/DNA/MegaMuga/")
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names = FALSE, colClasses="character")

# Remove individuals with > 20 % missing data
idx <- names(which(apply(genotypes,2,function(x){sum(is.na(x))}) /nrow(genotypes) > 0.2))
genotypes <- genotypes[,-which(colnames(genotypes) %in% idx)]

# Create numeric genotypes (not A, H, B)
genoNum <- apply(genotypes, 1, function(x){ as.numeric(as.factor(x)) })
rownames(genoNum) <- colnames(genotypes)
colnames(genoNum) <- rownames(genotypes)

corG <- cor(t(genoNum), use="pair")       # Calculate the correlation matrix

write.table(corG, "Analysis/genotypecorrelation.txt", sep = "\t", quote = FALSE)

pcas <- principal(corG, nfactors = 1)     # Calculate the first principal component

plot(pcas$loadings[,1])

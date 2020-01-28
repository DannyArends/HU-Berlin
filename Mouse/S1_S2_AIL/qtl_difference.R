# Modelling using AIL BFMI S1xS2
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written december, 2019

setwd("D:/Edrive/Mouse/S1_S2")

genotypes <- read.csv("genotypes.cleaned.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes_Jan2020.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
annotation <- read.csv("map.cleaned.txt", header=TRUE, sep="\t", check.names=FALSE)
annotation <- annotation[, c(1,2,3,6)]
colnames(annotation) <- c("Chromosome", "Position", "GenTrain Score", "SNP")
colnames(genotypes) <- gsub("V 888-", "" , colnames(genotypes)) 
rownames(phenotypes) <- paste0("AIL", rownames(phenotypes))

#Phenotype outliers get rid of them
noPhe <- c("Sex", "ID", "WG","Mutter","Grandma")

phenames <- colnames(phenotypes)[which(!(colnames(phenotypes) %in% noPhe))]

outliers <- apply(phenotypes[, phenames],2,function(x){
  up <- mean(x, na.rm=TRUE) + 3 * sd(x, na.rm=TRUE)
  low <- mean(x, na.rm=TRUE) - 3 * sd(x, na.rm=TRUE)
  x < low | x > up
})

for(x in phenames){
  idx <- which(outliers[,x])
  if(length(idx) > 0) phenotypes[idx,x] <- NA
}

# Scale down the phenotypes to the genotypes individuals
phenotypes <- phenotypes[colnames(genotypes),]

# Covariates we could/need to include in the model, we test them on their pvalue
wg <- as.numeric(phenotypes[, "WG"] == 8)

for(p in 57){
  mmatrix <- c()
  for(m in 1:nrow(genotypes)){
    a <- (as.numeric(factor(genotypes[m,], levels = c("A", "H", "B"))) - 2)
    d <- (as.numeric(genotypes[1,] == "H"))
    f <- factor(genotypes[m,], levels = c("A", "H", "B"))
    
    mdata <- data.frame(y = phenotypes[, p], WG = as.factor(wg), A = a, D = d, F = as.factor(f))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if(length(isNA) > 0) mdata <- mdata[-isNA, ]
    
    m2 <- lm(y ~ F, data = mdata)
    m1 <- lm(y ~ A + D, data = mdata)
    m0 <- lm(y ~ 1, data = mdata)
    pOA <- anova(m1, m0)
    pAD <- anova(m1)
    pF <- anova(m2, m0)

    effects <- coefficients(m1)[c("A", "D")]
    line <- c(pOA[2, "Pr(>F)"], as.numeric(na.omit(pAD[, "Pr(>F)"]))[1:2], pF[2, "Pr(>F)"], effects)
    mmatrix <- rbind(mmatrix, line)
  }
  colnames(mmatrix) <- c("pOA", "pWG", "pADD", "pDomDev", "pF", "Eff(ADD)", "Eff(DomDev)")
  rownames(mmatrix) <- rownames(genotypes)
}

-log10(mmatrix["UNCHS019508",])

for(p in 3:ncol(phenotypes)){
  mmatrix <- c()
  for(m in 1:nrow(genotypes)){
    F <- factor(genotypes[m,], levels = c("A", "H", "B"))
    
    mdata <- data.frame(cbind(y = phenotypes[, p], wg = wg, F = F))
    isNA <- which(apply(apply(mdata,1,is.na),2,any))
    if(length(isNA) > 0) mdata <- mdata[-isNA, ]
    
    m1 <- lm(y ~ wg + F, data = mdata)
    m0 <- lm(y ~ wg, data = mdata)
    pOA <- anova(m1, m0)
    pAD <- anova(m1)

    line <- c(pOA[2, "Pr(>F)"], as.numeric(na.omit(pAD[, "Pr(>F)"]))[2])
    mmatrix <- rbind(mmatrix, line)
  }
  colnames(mmatrix) <- c("pOA", "pF")
  rownames(mmatrix) <- rownames(genotypes)
}

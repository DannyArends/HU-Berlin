
library(parallel)

emma.kinship <- function(snps, method="additive", use="all") {
  n0 <- sum(snps==0,na.rm=TRUE)
  nh <- sum(snps==0.5,na.rm=TRUE)
  n1 <- sum(snps==1,na.rm=TRUE)
  nNA <- sum(is.na(snps))

  stopifnot(n0+nh+n1+nNA == prod(dim((snps))))

  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1

  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}

setwd("~/PopStructCow")

genotypes <- read.table("~/NAS/Cattle/DSN/analysis/QTL/genotypes.txt", sep="\t", check.names=FALSE)
phenotypes <- read.table("~/NAS/Cattle/DSN/analysis/QTL/phenotypes.txt", sep="\t")

DSNind <- rownames(phenotypes)[which(phenotypes[,"breed"] == "DSN")]
HFind <- rownames(phenotypes)[which(phenotypes[,"breed"] == "Holstein")]

genotypes.DSN <- genotypes[,DSNind]
genoE.DSN <- genotypes.DSN
genoE.DSN[genoE.DSN ==  0] <- 0.5
genoE.DSN[genoE.DSN == -1] <- 0

IBS.DSN <- emma.kinship(genoE.DSN)
colnames(IBS.DSN) <- rownames(IBS.DSN) <- colnames(genoE.DSN)
write.table(IBS.DSN, "IBS.DSN.txt", sep="\t")

genotypes.HF <- genotypes[,HFind]
genoE.HF <- genotypes.HF
genoE.HF[genoE.HF ==  0] <- 0.5
genoE.HF[genoE.HF == -1] <- 0

IBS.HF <- emma.kinship(genoE.HF)
colnames(IBS.HF) <- rownames(IBS.HF) <- colnames(genoE.HF)
write.table(IBS.HF, "IBS.HF.txt", sep="\t")

## No adjustments

cl <- makeCluster(getOption("cl.cores", 20))
pvals.HF <- unlist(parLapply(cl, rownames(genotypes), function(mname, genotypes, phenotypes, phe = "RZM"){
  geno  <- unlist(genotypes[mname,])
  pheno <- phenotypes[colnames(genotypes), phe]
  mymodel <- lm(pheno ~ geno)
  myanova <- anova(mymodel)
  return(myanova["Pr(>F)"]["geno",])
}, genotypes = genotypes.HF, phenotypes = phenotypes, phe = "RZM"))
pvals.DSN <- unlist(parLapply(cl, rownames(genotypes), function(mname, genotypes, phenotypes, phe = "RZM"){
  geno  <- unlist(genotypes[mname,])
  pheno <- phenotypes[colnames(genotypes), phe]
  mymodel <- lm(pheno ~ geno)
  myanova <- anova(mymodel)
  return(myanova["Pr(>F)"]["geno",])
}, genotypes = genotypes.DSN, phenotypes = phenotypes, phe = "RZM"))
stopCluster(cl)

lambda.HF <- median(qchisq(1.0 - pvals.HF, 1)) /  qchisq(0.5, 1)
lambda.DSN <- median(qchisq(1.0 - pvals.DSN, 1)) /  qchisq(0.5, 1)
cat("λ =", lambda.HF, lambda.DSN, "\n")

## Adjustment using the IBS matrix

IBS.HF <- 1-IBS.HF
HFcluster <- hclust(as.dist(IBS.HF))
fit.HF <- NULL
for(k in 2:(ncol(IBS.HF) / 5)){
  grouping <- cutree(HFcluster, k)
  samples <- names(grouping)
  grouping <- as.factor(grouping)
  myanova <- anova(lm(phenotypes[samples, "RZM"] ~ grouping))
  sqExplained <-  myanova["Sum Sq"]["grouping",]
  sqTotal <-  sum(myanova["Sum Sq"])
  fit.HF <- rbind(fit.HF, c(k, myanova["Pr(>F)"]["grouping",], sqExplained / sqTotal))
}
HF.top <- fit.HF[which.min(fit.HF[,2]), 1]

IBS.DSN <- 1-IBS.DSN
DSNcluster <- hclust(as.dist(IBS.DSN))
fit.DSN <- NULL
for(k in 2:(ncol(IBS.DSN) / 5)){
  grouping <- cutree(DSNcluster, k)
  samples <- names(grouping)
  grouping <- as.factor(grouping)
  myanova <- anova(lm(phenotypes[samples, "RZM"] ~ grouping))
  sqExplained <-  myanova["Sum Sq"]["grouping",]
  sqTotal <-  sum(myanova["Sum Sq"])
  fit.DSN <- rbind(fit.DSN, c(k, myanova["Pr(>F)"]["grouping",], sqExplained / sqTotal))
}
DSN.top <- fit.DSN[which.min(fit.DSN[,2]), 1]

cl <- makeCluster(getOption("cl.cores", 20))
pvals.HF.adj.nRO <- unlist(parLapply(cl, rownames(genotypes), function(mname, genotypes, phenotypes, subpop, phe = "ZWFPR"){
  geno  <- unlist(genotypes[mname,])
  pheno <- phenotypes[colnames(genotypes), phe]
  #mymodel <- lm(pheno ~ as.factor(subpop) + geno)
  mymodel <- lm(terms(pheno ~ as.factor(subpop) + geno, keep.order = TRUE))
  myanova <- anova(mymodel)
  return(myanova["Pr(>F)"]["geno",])
}, genotypes = genotypes.HF, phenotypes = phenotypes, subpop = cutree(HFcluster, HF.top), phe = "ZWFPR"))
pvals.DSN.adj <- unlist(parLapply(cl, rownames(genotypes), function(mname, genotypes, phenotypes, subpop, phe = "ZWFPR"){
  geno  <- unlist(genotypes[mname,])
  pheno <- phenotypes[colnames(genotypes), phe]
  mymodel <- lm(pheno ~ as.factor(subpop) + geno)
  myanova <- anova(mymodel)
  return(myanova["Pr(>F)"]["geno",])
}, genotypes = genotypes.DSN, phenotypes = phenotypes, subpop = cutree(DSNcluster, DSN.top), phe = "ZWFPR"))
stopCluster(cl)

lambda.HF <- median(qchisq(1.0 - pvals.HF.adj, 1)) /  qchisq(0.5, 1)
lambda.DSN <- median(qchisq(1.0 - pvals.DSN.adj, 1)) /  qchisq(0.5, 1)
cat("λ =", lambda.HF, lambda.DSN, "\n")


cl <- makeCluster(getOption("cl.cores", 20))
pheno <- phenotypes[colnames( genotypes.HF), phe]
subpop = as.factor(cutree(HFcluster, HF.top))
myadj <- lm(pheno ~ subpop)
npheno <- rep(NA, length(pheno))
nvalues <- myadj$residuals + mean(pheno)
npheno[as.numeric(names(nvalues))] <- nvalues

pvals.HF.all <- parLapply(cl, rownames(genotypes), function(mname, genotypes, pheno, phe = "ZWFPR"){
  geno  <- unlist(genotypes[mname,])
  mymodel <- lm(pheno ~ geno)
  myanova <- anova(mymodel)
  return(c(myanova["Pr(>F)"]["subpopf",], myanova["Pr(>F)"]["geno",], myanova["Pr(>F)"]["subpopf:geno",]))
}, genotypes = genotypes.HF, pheno = npheno, phe = "ZWFPR")
stopCluster(cl)


mTop <- 17995
phe = "ZWFPR"
subpop = cutree(HFcluster, HF.top)
 
geno  <- unlist( genotypes.HF[mTop,])
pheno <- phenotypes[colnames( genotypes.HF), phe]
subpopf <- as.factor(subpop)
mymodel <- lm(terms(pheno ~ subpopf + geno, keep.order = TRUE))
myanova.top <- anova(mymodel)

mRand <- 17

geno  <- unlist( genotypes.HF[mRand,])
pheno <- phenotypes[colnames( genotypes.HF), phe]
subpopf <- as.factor(subpop)
mymodel <- lm(terms(pheno ~ subpopf + geno, keep.order = TRUE))
myanova.rnd <- anova(mymodel)




geno.rnd <- unlist( genotypes.HF[mRand,])
geno.top <- unlist( genotypes.HF[mTop,])

mymodel.rnd <- lm(terms(npheno ~ geno.rnd, keep.order = TRUE))
mymodel.top <- lm(terms(npheno ~ geno.top, keep.order = TRUE))
myanova.rnd <- anova(mymodel.rnd)
myanova.top <- anova(mymodel.top)



startChr14 <- which(grepl("Chr14", rownames(genotypes.HF)))[1:100]
corM <- cor(t(genotypes.HF[startChr14,]), use="pair", method="spear")
lodscores <- -log10(unlist(lapply(pvals.HF.all,"[",2)))
lodscoresB <- -log10(pvals.HF.adj.nRO)

topM <- rownames(genotypes.HF)[which.max(-log10(unlist(lapply(pvals.HF.all,"[",2))))]

plot(c(0,1), c(0,150),t='n')
points(abs(corM[topM,]), lodscores[startChr14])
points(abs(corM[topM,]), lodscoresB[startChr14], col="orange")
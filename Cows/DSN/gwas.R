library(parallel)

setwd("D:/Edrive/Cow/DSN")

genotypes <- read.table("genotypes.txt", sep="\t", check.names=FALSE)
phenotypes <- read.table("phenotypes.txt", sep="\t")

# order by genotypes
phenotypes <- phenotypes[colnames(genotypes),]
phenotype <- "MKG3_1"

genoCorM <- matrix(0, ncol(genotypes), ncol(genotypes))
for(x in 1:10){
  genoCorM <- genoCorM + cor(genotypes[sample(nrow(genotypes), 2500), ], use="pair")
  cat(x, "\n")
}
genoCorM <- (genoCorM / 10)

gts <- apply(genotypes, 1, table)
mingts <- lapply(gts, min)

genoclean <- t(apply(genotypes, 1, function(x){
  tbl <- table(x)
  if(min(tbl) < 30){
    toosmall <- names(tbl)[which.min(tbl)]
    x[x==toosmall] <- NA
  }
  return(x)
}))

yob <- as.factor(phenotypes[,"GEBJ"])
plot(phenotypes[,phenotype] ~ yob)
father <- as.factor(phenotypes[,"VATID"])
plot(phenotypes[,phenotype] ~ father)
farm <- as.factor(phenotypes[,"BETR"])
plot(phenotypes[,phenotype] ~ farm)

cl <- makeCluster(getOption("cl.cores", 6))
clusterExport(cl, "yob") # export yob to the nodes
clusterExport(cl, "father") # export father to the nodes
clusterExport(cl, "farm") # export farm to the nodes

mapping <- parApply(cl, genotypes, 1,function(marker, phenotype){
  addeff <- marker
  domdev <- as.numeric(marker == 1)
  mymodel <- lm(phenotype ~ yob + father + farm + addeff)
  myanova <- anova(mymodel)
  return(myanova)
}, phenotype = phenotypes[,phenotype])
stopCluster(cl)

values = lapply(mapping, function(x){return(x[[5]])})
pvals = matrix(unlist(values), length(mapping), length(mapping[[1]][[5]]), byrow=TRUE)

expected <- (rank(pvals[,4], ties.method="first")+0.5) / (length(pvals[,4]) + 1)
plot(-log10(expected), -log10(pvals[,4]))

lambda <- round(median(qchisq(1.0 - pvals[,4], 1), na.rm=TRUE) /  qchisq(0.5, 1),3)

cl <- makeCluster(getOption("cl.cores", 6))
clusterExport(cl, "yob") # export yob to the nodes
clusterExport(cl, "father") # export father to the nodes
clusterExport(cl, "farm") # export farm to the nodes

mapping2 <- parApply(cl, genoclean, 1, function(marker, phenotype){
  addeff <- marker
  domdev <- as.numeric(marker == 1)
  mymodel <- lm(phenotype ~ yob + father + farm + addeff)
  myanova <- anova(mymodel)
  return(myanova)
}, phenotype = phenotypes[,phenotype])
stopCluster(cl)

values2 = lapply(mapping2, function(x){return(x[[5]])})
pvals2 = matrix(unlist(values2), length(mapping2), length(mapping2[[1]][[5]]), byrow=TRUE)

expected2 <- (rank(pvals2[,4], ties.method="first")+0.5) / (length(pvals2[,4]) + 1)
plot(-log10(expected2), -log10(pvals2[,4]))

lambda2 <- round(median(qchisq(1.0 - pvals2[,4], 1), na.rm=TRUE) /  qchisq(0.5, 1),3)


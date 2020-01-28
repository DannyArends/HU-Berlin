
library(parallel)

setwd("D:/Edrive/Mouse/S1_S2")

genotypes <- read.table("genotypes.cleaned.txt", sep="\t", colClasses = "character")
map <- read.table("map.cleaned.txt", sep="\t")
phenotypes <- read.table("allPhenotypes.txt", sep="\t")

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, map[map[,"chr"] == chr,])
}

phenames <- colnames(phenotypes)[-c(1,2, 58,59)]

# Getting rid of the outliers
outliers <- apply(phenotypes[, phenames],2, function(x){
  up <- mean(x, na.rm = TRUE) + 3 * sd(x, na.rm = TRUE)
  low <- mean(x, na.rm=TRUE) - 3 * sd(x, na.rm=TRUE)
  x < low | x > up
 })
 
for(x in phenames){
  idx <- which(outliers[,x])
  if(length(idx) > 0) phenotypes[idx,x] <- NA
}

# Adjust the tissue weight for the total bodyweight
tissues <- colnames(phenotypes[, c(46:55)])
for (x in tissues){
  phenotypes[, x] <- phenotypes[, x] / phenotypes[, "Gewicht"]
}

# Make sure that the ordering between phenotypes and genotypes matches 
colnames(genotypes) <- gsub("AIL", "", colnames(genotypes))
phenotypes <- phenotypes[colnames(genotypes),]
genotypes <- genotypes[, rownames(phenotypes)]
write.table(genotypes, "OrderedGenotypes.txt", sep = "\t", quote=FALSE)

# QTL mapping
cl <- makeCluster(4)
pvals <- c()
for(x in "Triglycerides"){
  pvals <- rbind(pvals, 
    parApply(cl, genotypes, 1, function(gts, pheno){
      marker <- as.numeric(factor(as.character(gts), levels = c("A", "H", "B")))
      mylm <- lm(pheno ~ marker)                                                    # We don´t include mother or litter size as cofactors
      myanova <- anova(mylm)
      return(myanova[[5]][1])
    }, pheno = phenotypes[,x])
  )
  cat("Done", x, "\n")
}
stopCluster(cl)
rownames(pvals) <- phenames
lods <- -log10(pvals)
apply(lods,1,max)
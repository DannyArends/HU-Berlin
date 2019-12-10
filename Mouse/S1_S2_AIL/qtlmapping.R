# Basic QTL mapping
library(parallel)
setwd("D:/Edrive/Mouse/S1_S2")

genotypes <- read.table("genotypes.cleaned.txt", sep="\t", colClasses = "character")
map <- read.table("map.cleaned.txt", sep="\t")
phenotypes <- read.table("allPhenotypes.txt", sep="\t")

# Fix rownames
rownames(phenotypes) <- gsub("V 888-", "AIL", phenotypes[, "ID"])
phenames <- colnames(phenotypes)[3:(ncol(phenotypes)-2)]

outliers <- apply(phenotypes[, phenames],2,function(x){
  up <- mean(x, na.rm=TRUE) + 3 * sd(x, na.rm=TRUE)
  low <- mean(x, na.rm=TRUE) - 3 * sd(x, na.rm=TRUE)
  x < low | x > up
})

for(x in phenames){
  idx <- which(outliers[,x])
  if(length(idx) > 0) phenotypes[idx,x] <- NA
}

#Additive: -1,0,1
#Dom A: 0,0,1
#Dom B: 1,0,0
#DomDev: -1,0,1 + 0,1,0

# Take only animals with genotypes
phenotypes <- phenotypes[colnames(genotypes),]

cl <- makeCluster(6)
pvals <- c()
for(x in "Triglycerides"){
  pvals <- rbind(pvals, 
    parApply(cl, genotypes, 1, function(gts, pheno){
      marker <- as.numeric(factor(as.character(gts), levels = c("A", "H", "B")))
      mylm <- lm(pheno ~ marker)
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

library(RColorBrewer)
colz <- brewer.pal(n = 9, name = "PuRd")
chrs <- 1:21 
names(chrs) <- c(1:19, "X", "Y")

for(x in "Triglycerides"){
  plot(c(1,21), c(0,200000000), t = 'n', xaxt = "n", las= 2, ylab = "Position (mb)", xlab = "Chr", yaxt = 'n', main=x)
  cnt <- 1
  aa <- apply(map, 1, function(r) { 
    points(x = chrs[r[1]], y = r[2], pch = "-", col = colz[ceiling(lods[x, cnt])], cex=3); 
    cnt <<- cnt + 1
    }
  )
  axis(1, at = chrs, names(chrs))
  axis(2, at = seq(0,200000000, 25000000), seq(0,200000000, 25000000)/1000000)
  readline(prompt="Press [enter] to continue")
}

maxMarker <- function(phe){
  cat("Max marker:", names(which.max(lods[phe,])), "\n")
  return(factor(as.character(genotypes[names(which.max(lods[phe,])),]), levels = c("A", "H", "B")))
}



boxplot(phenotypes[, "Triglycerides"] ~ maxMarker("Triglycerides"))
boxplot(phenotypes[, "Leber"] ~ factor(as.character(genotypes[names(which.max(lods["Leber",])),]), levels = c("A", "H", "B")))

#model trygl
m1 <- factor(genotypes["JAX00632487",], levels = c("A", "H", "B"))
m2 <- factor(genotypes["UNC27577908",], levels = c("A", "H", "B"))
m3 <- factor(genotypes["JAX00064248",], levels = c("A", "H", "B"))

anova(lm(phenotypes[, "Triglycerides"] ~ m1 + m2 * m3))

#Day 126

maxMarker("D126")

m1 <- factor(genotypes["UNCHS040893",], levels = c("A", "H", "B"))
m2 <- factor(genotypes["UNC415288",], levels = c("A", "H", "B"))

anova(lm(phenotypes[, "D126"] ~ m1 + m2))

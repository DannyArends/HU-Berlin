#
# Analyze all horse data
#

setwd("E:/Horse/DNA/")

genotypes_snp  <- read.csv("combined/input/genotypes_snp.txt", sep="\t")
markerinfo     <- read.csv("combined/input/map.txt", sep="\t")
genotypes_num  <- read.csv("combined/input/genotypes_num.txt", sep = "\t")
genotypes_AB   <- read.csv("combined/input/genotypes_AB.txt", sep = "\t")
phenotypes     <- read.csv("combined/input/phenotypes.txt", sep="\t")

# How many markers on each chromosome
for(x in unique(markerinfo[,"Chr"])){
  cat(x, length(which(markerinfo[,"Chr"] == x)),"\n")
}

ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K","H","S", "PrzH", "Arab", "Arabian", "Norwegian Fjord"))

strains <- as.character(phenotypes[ii,"Strain"])
names(strains) <- rownames(phenotypes)[ii]

# Create colors
cols <- c("red", "blue", "orange", "black", "brown", "brown", "purple")
names(cols) <- c("K","H","S", "PrzH", "Arab", "Arabian", "Norwegian Fjord")

labelCol <- function(x) {
  if (is.leaf(x)) {
    hclass <- strains[attr(x, "label")]             # Fetch the class label
    hcol <- cols[hclass]                            # Determine color of the label
    cat(attr(x, "label"), hclass, hcol, "\n")
    attr(x, "nodePar") <- list(lab.col=hcol)
  }
  return(x)
}

clusters <- hclust(dist(t(genotypes_num[,rownames(phenotypes)[ii]]),"manhattan"))
dendrogram <- as.dendrogram(clusters)
dendrogram.col <- dendrapply(dendrogram, labelCol)
plot(dendrogram.col, main = "")

library(StAMPP)

ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K", "H", "S", "PrzH", "Norwegian Fjord"))
strains <- as.character(phenotypes[ii,"Strain"])
stammpinput <- t(genotypes_AB[, rownames(phenotypes)[ii]])

stammpinput <- data.frame(cbind(rownames(stammpinput), strains, 2, "BiA", stammpinput))
colnames(stammpinput)[1:4] <- c("Sample", "Pop", "Ploidy", "Format")

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values

stammpinput.fst <- stamppFst(stammpinput.freq, 1000, 95, 4) # Population Fst values

## GWAS

ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K","H","S"))

# Correlation matrix over classical phenotypes
cor(phenotypes[,c("WH", "CW", "CH", "NG", "TG", "ChG", "ChD", "ChW", "BLL", "BL", "FCL", "HCL", "Distance..km.", "Speed.km.hr.")], use="pair", method="spearman")

# Analyse the different covariates for the different phenotypes
sex     <- as.factor(unlist(phenotypes[ii,"Sex"]))
strain  <- as.factor(unlist(phenotypes[ii,"Strain"]))
ageAtM  <- as.numeric(unlist(phenotypes[ii,"D.Measure"])) - as.numeric(unlist(phenotypes[ii,"D.Birth"]))
ageAtR  <- as.numeric(unlist(phenotypes[ii,"D.Racing"])) - as.numeric(unlist(phenotypes[ii,"D.Birth"]))

classicalpheno <- c("WH", "CW", "CH", "NG", "TG", "ChG", "ChD", "ChW", "BLL", "BL", "FCL", "HCL")

racepheno <- c("Distance..km.", "Speed.km.hr.")

# Calculate thesignificance of the fixed effects
pClassical <- matrix(NA, length(classicalpheno), 3, dimnames = list(classicalpheno, c("Sex", "Age", "Strain")))
for(ph in classicalpheno) {
  model <- anova(lm(as.numeric(phenotypes[ii,ph]) ~ sex + ageAtM + strain))
  pClassical[ph,] <- model[[5]][1:3]
}

pRace <- matrix(NA, length(racepheno), 3, dimnames=list(racepheno, c("Sex", "Age", "Strain")))
for(ph in racepheno) {
  model <- anova(lm(as.numeric(phenotypes[ii,ph]) ~ sex + ageAtR + strain))
  pRace[ph,] <- model[[5]][1:3]
}

write.table(rbind(pClassical, pRace), "combined/output/covariates.txt", sep="\t")

# Calculate the phenotypes after adjusting the fixed effects
phenoC <- matrix(NA, length(c(classicalpheno, racepheno)), nrow(phenotypes[ii,]), dimnames = list(c(classicalpheno, racepheno), rownames(phenotypes[ii,])))
for(ph in classicalpheno) {
  model <- lm(as.numeric(phenotypes[ii,ph]) ~ sex + ageAtM + strain)
  pCorrected <- model$residuals + model$coefficients["(Intercept)"]
  phenoC[ph, as.numeric(names(pCorrected))] <- pCorrected
}
for(ph in racepheno) {
  model <- lm(as.numeric(phenotypes[ii,ph]) ~ sex + ageAtR + strain)
  pCorrected <- model$residuals + model$coefficients["(Intercept)"]
  phenoC[ph, as.numeric(names(pCorrected))] <- pCorrected
}

genotypes <- genotypes_num[,ii]

# Do the GWAS using a single QTL model on the data adjusted for fixed effects
pvalues <- matrix(NA, length(c(classicalpheno, racepheno)), nrow(genotypes), dimnames = list(c(classicalpheno, racepheno), rownames(genotypes)))
for(phe in rownames(phenoC)) {
  cat("Computing GWAS results for:", phe, "\n")
  pvalues[phe, ] <- apply(genotypes, 1, function(marker) {
    tryCatch(res <- anova(lm(as.numeric(phenoC[phe,]) ~ as.factor(marker)))[[5]][1], error = function(e){ res <<- NA })
    return(res)
  })
}
pvalues <- t(pvalues)
write.table(pvalues, "combined/output/pvaluesGWAS.txt", sep="\t")

pvalues     <- read.table("combined/output/pvaluesGWAS.txt", sep="\t")
markerinfo  <- read.csv("combined/input/map.txt", sep="\t")

# Summarize significant results in a file
results <- NULL
for(phe in colnames(pvalues)) {
  ii <- which(pvalues[, phe] < (0.05/nrow(pvalues)))
  if(length(ii) > 0){
    for(i in ii){
      results <- rbind(results, c(phe, rownames(pvalues)[i], as.character(markerinfo[as.character(rownames(pvalues)[i]),"Chr"]), markerinfo[as.character(rownames(pvalues)[i]),"MapInfo"], pvalues[i, phe]))
    }
    cat(phe, rownames(pvalues)[ii],"\n")
  }
}
write.table(results, "combined/output/pvaluesGWAS_0.05.txt", sep="\t", col.names=FALSE, row.names=FALSE)


pvalues     <- read.table("combined/output/pvaluesGWAS.txt", sep="\t")
markerinfo  <- read.csv("combined/input/map.txt", sep="\t")


# Create the manhattan plots
neworder <- order(markerinfo[,"MapInfo"])
markerinfo <- markerinfo[neworder, ]
pvalues <- pvalues[neworder, ]

newmarkerinfo <- NULL
newpvalues <- NULL

for(x in c("1", "2", "3", "4", "5", "6", "7", "8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "X")){
  msubset <- markerinfo[which(markerinfo[, "Chr"] == x),]
  newmarkerinfo <- rbind(newmarkerinfo, msubset)
  newpvalues <- rbind(newpvalues, pvalues[rownames(msubset),])
}

markerinfo <- newmarkerinfo
pvalues <- newpvalues

chrcols <- 1+ (as.numeric(as.factor(markerinfo[,"Chr"])) %% 2)

for(phe in colnames(pvalues)) {
  png(paste0("combined/output/GWAS_", gsub("\\\\hr", "pH", phe), ".png"), width=1024, height=600)
    plot(x = c(1, nrow(pvalues)), y = c(0,10), t='n', main=phe, ylab = "LOD score (-log10(p))", xlab = "Marker")
    points(-log10(pvalues[,phe]), pch = 19, cex = 0.5, col=chrcols)
    abline(h = -log10(0.1/nrow(pvalues)), col="orange")
    abline(h = -log10(0.05/nrow(pvalues)), col="gold")
    abline(h = -log10(0.01/nrow(pvalues)), col="green")
  dev.off()
}

#Principal component analysis
library(pcaMethods)
resSvd <- pca(genotypes_num, method = "svd", nPcs = 5, center = FALSE)

colz <- as.numeric(as.factor(phenotypes[rownames(resSvd@loadings),"Strain"])) + 1
cc <- as.numeric(as.factor(unique(phenotypes[rownames(resSvd@loadings),"Strain"]))) + 1
names(cc) <- unique(phenotypes[rownames(resSvd@loadings),"Strain"])

plot(resSvd@loadings[,1:2], col=colz, pch=18)

table(phenotypes[which(!phenotypes[,"Strain"] %in% c("K","H","S","Kab", "Arab", "EV", "PrzH", "ExP", "AnK")), "Strain"]
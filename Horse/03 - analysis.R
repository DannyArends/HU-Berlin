#
# Analyze all horse data
#

toGenPop <- function(genotypes){
  numericG <- apply(genotypes, 1, function(x){
    geno <- table(unlist(strsplit(as.character(x),"")))
    #cat(names(geno),"\n")
    a1 <- paste0(names(geno)[1],names(geno)[1])
    a2 <- paste0(sort(c(names(geno)[1],names(geno)[2])),collapse="")
    a3 <- paste0(names(geno)[2],names(geno)[2])
    ngeno <- rep(NA,length(x))
    ngeno[x == a1] <- "0101"
    ngeno[x == a2] <- "0102"
    ngeno[x == a3] <- "0202"
    return(ngeno)
  })
  rownames(numericG) <- colnames(genotypes)
  return(t(numericG))
}

setwd("D:/Edrive/Horse/DNA/")

genotypes_snp  <- read.csv("combined/input/genotypes_snp.txt", sep="\t")
markerinfo     <- read.csv("combined/input/map.txt", sep="\t")
genotypes_num  <- read.csv("combined/input/genotypes_num.txt", sep = "\t")
genotypes_AB   <- read.csv("combined/input/genotypes_AB.txt", sep = "\t")
phenotypes     <- read.csv("combined/input/phenotypes.txt", sep="\t")

genotypes_genpop <- t(toGenPop(genotypes_AB))
set.seed(0)
#rsample <- sample(ncol(genotypes_genpop), ncol(genotypes_genpop) / 20)
#genotypes_genpop <- genotypes_genpop[, rsample]

# How many markers on each chromosome
#for(x in unique(markerinfo[,"Chr"])){
#  cat(x, length(which(markerinfo[,"Chr"] == x)),"\n")
#}

ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K","H","S", "PrzH", "Arab", "Arabian", "Norwegian Fjord"))

strains <- as.character(phenotypes[ii,"Strain"])
names(strains) <- rownames(phenotypes)[ii]

## Write the genepop stucture to disk

genotypes_genpop <- genotypes_genpop[which(rownames(genotypes_genpop) %in% names(strains)),]
rownames(genotypes_genpop) <- gsub(" ","", paste0(strains[rownames(genotypes_genpop)],","))

cat("BLANK\n", file="combined/input/genotypes_genpop.txt")
for(x in colnames(genotypes_genpop)){
  cat(paste0(x, "\n"), file="combined/input/genotypes_genpop.txt", append=TRUE)
}
for(pop in unique(rownames(genotypes_genpop))) {
  cat("POP\n", file="combined/input/genotypes_genpop.txt", append=TRUE)
  ii <- which(rownames(genotypes_genpop) == pop)
  write.table(genotypes_genpop[ii,], sep = " ", quote=FALSE, na = "0000", file="combined/input/genotypes_genpop.txt", append=TRUE, col.names=FALSE)
}
cat("POP\n", file="combined/input/genotypes_genpop.txt", append=TRUE)
ii <- which(rownames(genotypes_genpop) %in% c("S,", "K,", "H,"))
genotypes_copy <- genotypes_genpop[ii,]
rownames(genotypes_copy) <- rep("CombinedSHK", nrow(genotypes_copy))
write.table(genotypes_copy[,], sep = " ", quote=FALSE, na = "0000", file="combined/input/genotypes_genpop.txt", append=TRUE, col.names=FALSE)


ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K","H","S", "PrzH"))

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

usingsamples <- rownames(phenotypes)[ii]
usingsamples <- usingsamples[-which(usingsamples == "ARR1")]


clusters <- hclust(dist(t(genotypes_num[,usingsamples]),"manhattan"))
dendrogram <- as.dendrogram(clusters)
dendrogram.col <- dendrapply(dendrogram, labelCol)
plot(dendrogram.col, main = "Distance (manhattan)")


ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K","H","S", "PrzH", "Arabian"))

strains <- as.character(phenotypes[ii,"Strain"])
names(strains) <- rownames(phenotypes)[ii]

# Create colors
cols <- c("red", "blue", "orange", "black", "brown", "purple")
names(cols) <- c("K","H","S", "PrzH", "Arabian", "Norwegian Fjord")
cnt <- 1
labelCol <- function(x) {
  if (is.leaf(x)) {
    hclass <- strains[attr(x, "label")]             # Fetch the class label
    hcol <- cols[hclass]                            # Determine color of the label
    cat(attr(x, "label"), hclass, hcol, "\n")
   # attr(x, "nodePar") <- list(col=hcol, cex=c(2,4))
    #attr(x, "leaflab") <- "none"
    if(grepl("NORF", attr(x, "label"))){
      #attr(x, "label") <- paste0("NORF",cnt)
      cnt <<- cnt + 1
    }
    attr(x, "label") <- strsplit(hclass, "")[[1]][1]
  }
  return(x)
}

#clusters <- hclust(dist(t(genotypes_num[,rownames(phenotypes)[ii]]),"manhattan"))
dendrogram <- as.dendrogram(clusters)
dendrogram.col <- dendrapply(dendrogram, labelCol)
postscript("dendrogramArabPrzNorf_29_8.eps", width = 16.0, height = 4.0, horizontal = FALSE, onefile = FALSE, paper = "special")
plot(dendrogram.col, main = "", las=2, leaflab ="none", horiz = FALSE)
text(1:length(labels(dendrogram.col)), -400, labels = labels(dendrogram.col))
dev.off()

### stampp ANALYSIS
library(StAMPP)

ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K", "H", "S", "PrzH", "Norwegian Fjord", "Arabian"))
strains <- as.character(phenotypes[ii,"Strain"])
stammpinput <- t(genotypes_AB[, rownames(phenotypes)[ii]])

stammpinput <- data.frame(cbind(rownames(stammpinput), strains, 2, "BiA", stammpinput))
colnames(stammpinput)[1:4] <- c("Sample", "Pop", "Ploidy", "Format")

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values
stammp.D.pop1 <- stamppNeisD(stammpinput.freq, FALSE) # Population D values

stammpamova1 <- stamppAmova(stammp.D.pop1, stammpinput.freq)

stammpinput.fst <- stamppFst(stammpinput.freq, 1000, 95, 4) # Population Fst values

### diveRsity ANALYSIS

install.packages("diveRsity")

setwd("E:/Horse/DNA/")
library(diveRsity)
starttime <-  proc.time()

basicStats <- divBasic(infile = "combined/input/genotypes_genpop.txt", outfile="fstOnlyOut.txt", gp=2, bootstraps = 100)
names(basicStats$fis) <- colnames(basicStats$Ho)

endtime <-  proc.time()

basicStats$Ho["overall",] # Observed
basicStats$He["overall",] # Expected
basicStats$fis[["S,"]]["overall",]              # Fis
basicStats$fis[["H,"]]["overall",]              # Fis
basicStats$fis[["K,"]]["overall",]              # Fis
basicStats$fis[["CombinedSHK"]]["overall",]     # Fis

advancedStats <- diffCalc(infile = "combined/input/genotypes_genpop.txt", outfile="fstOnlyOut.txt", fst=TRUE, pairwise=TRUE, boots = 1000)
# Pairwise Fst per populations
advancedStats$pairwise$Fst


saklawi <- names(which(strains == "S"))
aaa <- apply(genotypes_AB[,saklawi], 1, function(x){sum(x == "AB",na.rm=TRUE) / sum(!is.na(x))})
inBS <- names(basicStats$Ho[,"S,"])
plot(aaa[inBS], basicStats$Ho[,"S,"])


genotypes_AB["BIEC2_506494", saklawi]
aaa["BIEC2_506494"]
basicStats$Ho["BIEC2_506494","S,"]
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

setwd("D:/Edrive/Horse/DNA/")
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



mlenghts <- table(unlist(as.numeric(as.character(newmarkerinfo[,"Chr"]))))
mlenghtsd2 <- table(unlist(as.numeric(as.character(newmarkerinfo[,"Chr"])))) / 2

chrend <- c()
chrstart <- c(0)
for(x in 1:32){
  chrend <- c(chrend, chrstart[x] + mlenghts[x])
  chrstart <- c(chrstart, chrend[x])
}
names(chrstart) <- 1:32
names(chrend) <- 1:32

markerinfo <- newmarkerinfo
pvalues <- newpvalues

chrcols <- 1 + (as.numeric(as.character(markerinfo[,"Chr"])) %% 2)
chrcols[is.na(chrcols)] <- 1

chrcolsA <- as.numeric(as.character(markerinfo[,"Chr"]))
chrcolsA[is.na(chrcolsA)] <- 1

for(phe in colnames(pvalues)) {
  png(paste0("combined/output/GWAS_", gsub("\\\\hr", "pH", phe), ".png"), width=1024, height=600)
   # op <- par(mfrow = c(2,1))
   # plot(x = c(1, nrow(pvalues)), y = c(0,10), t='n', main=phe, ylab = "LOD score (-log10(p))", xlab = "Marker")
   # points(-log10(pvalues[,phe]), pch = 19, cex = 0.5, col=chrcols)
   #  abline(v=chrstart)
     
    plot(x = c(1, nrow(pvalues)), y = c(0,10), t='n', main=phe, ylab = "LOD score (-log10(p))", xlab = "Marker")
    points(-log10(pvalues[,phe]), pch = 19, cex = 0.5, col=chrcolsA)
    #abline(v=chrstart)
    abline(h = -log10(0.1/nrow(pvalues)), col="orange")
    abline(h = -log10(0.05/nrow(pvalues)), col="gold")
    abline(h = -log10(0.01/nrow(pvalues)), col="green")
  dev.off()
}

#Principal component analysis

ii <- which(as.character(phenotypes[,"Strain"]) %in% c("K","H","S", "Arabian"))

strains <- as.character(phenotypes[ii,"Strain"])
names(strains) <- rownames(phenotypes)[ii]



geno_num <- genotypes_num[,which(colnames(genotypes_num) %in% names(strains))]

misdata <- which(apply(geno_num, 1, function(x){any(is.na(x))}))
geno_num <- geno_num[-misdata,]

res <- apply(t(geno_num),2,table)
markersBad <- names(which(lapply(res, length) == 1))
geno_num <- geno_num[-which(rownames(geno_num) %in% markersBad),]

pcares <- prcomp(t(geno_num), scale=TRUE)
sumpca <- summary(pcares)

groups <- strains

plot(pcares$x[,1], pcares$x[,2], col=as.numeric(as.factor(groups)), cex=1.6, pch=18)
legend("topleft", col=as.numeric(unique(as.factor(groups))), legend = unique(as.factor(groups)),cex = 0.8, pch=18)

# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}

# Variable correlation/coordinates
var.coord <- t(apply(pcares$rotation, 1, var_cor_func, pcares$sdev))
head(var.coord[, 1:4])
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))

highcontrib <- names(which(var.contrib[,1] > (100 / length(var.contrib[,1])) * 10))

markerinfo[highcontrib,]


mafs <- NULL
for(x in unique(groups)){
	inG <- names(groups)[which(groups == x)]
	geno_num[highcontrib,inG]
	mafs <- rbind(mafs, apply(geno_num[highcontrib,inG],1,function(x){
	  ss <- (sum(as.numeric(x == 1)) * 2) + sum(as.numeric(x == 2))
	  cat(ss, (2*sum(table(x))), "\n")
	  maf <- (ss / (2*sum(table(x)))) * 100
	  return(maf)
	}))
}
rownames(mafs) <- unique(groups)

cbind(markerinfo[highcontrib,c(10,11)], t(mafs))


library(pcaMethods)


resSvd <- pca(genotypes_num, method = "svd", nPcs = 5, center = FALSE)

colz <- as.numeric(as.factor(phenotypes[rownames(resSvd@loadings),"Strain"])) + 1
cc <- as.numeric(as.factor(unique(phenotypes[rownames(resSvd@loadings),"Strain"]))) + 1
names(cc) <- unique(phenotypes[rownames(resSvd@loadings),"Strain"])

plot(resSvd@loadings[,1:2], col=colz, pch=18)

table(phenotypes[which(!phenotypes[,"Strain"] %in% c("K","H","S","Kab", "Arab", "EV", "PrzH", "ExP", "AnK")), "Strain"]
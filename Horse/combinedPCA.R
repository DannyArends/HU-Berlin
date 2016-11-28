# Combined PCA analysis cPCA for population genetics for all the horse data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2016
# first written Feb, 2015

setwd("D:/Edrive/Horse/DNA/")

markerinfo <- read.csv("combined/input/map.txt", sep="\t")
genotypes_num <- read.csv("combined/input/genotypes_num.txt", sep = "\t")
phenotypes <- read.csv("combined/input/phenotypes.txt", sep="\t", colClasses="character")

phenotypes[phenotypes[,"Strain"] == "S", "Strain"] <- "Arabian"
phenotypes[phenotypes[,"Strain"] == "K", "Strain"] <- "Arabian"
phenotypes[phenotypes[,"Strain"] == "H", "Strain"] <- "Arabian"
phenotypes[phenotypes[,"Strain"] == "Kab", "Strain"] <- "Kabarda"

npStrain <- table(phenotypes[,"Strain"])
strains <- names(which(npStrain > 10))

# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){ return(var.loadings*comp.sdev) }
contrib <- function(var.cos2, comp.cos2){ return((var.cos2 * 100)/comp.cos2) }

misdata <- which(apply(genotypes_num, 1, function(x){any(is.na(x))}))
genotypes_num <- genotypes_num[-misdata, ]

M <- NULL
for (x in 1:length(strains)) {
  for (y in 1:length(strains)) {
    cat(strains[x], "versus", strains[y], "\n")
    indX <- rownames(phenotypes)[which(phenotypes[,"Strain"] == strains[x])]
    indY <- rownames(phenotypes)[which(phenotypes[,"Strain"] == strains[y])]
    geno_num <- genotypes_num[, c(indX,indY)]
    
    # Remove markers with missing data
    misdata <- which(apply(geno_num, 1, function(x){any(is.na(x))}))
    geno_num <- geno_num[-misdata, ]
    
    # Remove markers with a single genotype
    res <- apply(t(geno_num), 2, table)
    markersBad <- names(which(lapply(res, length) == 1))
    geno_num <- geno_num[-which(rownames(geno_num) %in% markersBad), ]
    
    cat("Using", nrow(geno_num), "markers\n")
    
    pcares <- prcomp(t(geno_num), scale=TRUE)
    sumpca <- summary(pcares)

    # Variable correlation/coordinates
    var.coord <- t(apply(pcares$rotation, 1, var_cor_func, pcares$sdev))
    var.cos2 <- var.coord^2
    comp.cos2 <- apply(var.cos2, 2, sum)
    var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
    top100 <- names(sort(var.contrib[,1], decreasing=TRUE)[1:100])
    M <- rbind(M, c(strains[x], strains[y], top100))
  }
}

write.table(M, "top100PCA.txt", sep = "\t")


setwd("D:/Edrive/Horse/DNA/")
markerinfo     <- read.csv("combined/input/map.txt", sep="\t")
genotypes_num  <- read.csv("combined/input/genotypes_num.txt", sep = "\t")
phenotypes     <- read.csv("combined/input/phenotypes.txt", sep="\t", colClasses="character")

phenotypes[phenotypes[,"Strain"] == "S", "Strain"] <- "Arabian"
phenotypes[phenotypes[,"Strain"] == "K", "Strain"] <- "Arabian"
phenotypes[phenotypes[,"Strain"] == "H", "Strain"] <- "Arabian"
phenotypes[phenotypes[,"Strain"] == "Kab", "Strain"] <- "Kabarda"

npStrain <- table(phenotypes[,"Strain"])
strains <- names(which(npStrain > 10))

M <- read.table("top100PCA.txt", sep = "\t", colClasses="character")

markerinfo <- markerinfo[with(markerinfo, order(Chr, MapInfo)), c("Chr", "MapInfo")]
dim(genotypes_num)
genotypes_num <- genotypes_num[rownames(markerinfo),]
dim(genotypes_num)

cols <- 1:32
names(cols) <- strains

chrs <- c(1:31,"X")

chrs.starts <- c("1" = 0)
chrs.lengths <- c()
for(x in chrs){
  chrs.lengths <- c(chrs.lengths, max(as.numeric(markerinfo[markerinfo[,"Chr"] == x, "MapInfo"])) + 25000000)
  chrs.starts <- c(chrs.starts, (chrs.starts[length(chrs.starts)] + chrs.lengths[length(chrs.lengths)]))
}
names(chrs.starts) <- chrs
names(chrs.lengths) <- chrs

res <- matrix(0, length(strains), nrow(genotypes_num), dimnames=list(strains, rownames(genotypes_num)))
for (x in 1:nrow(M)) {
  if(M[x,1] != M[x,2]){
    mNames <- as.character(M[x, 3:102])
    res[M[x,1], mNames] <- res[M[x,1], mNames] + 1
  }
}

#res <- res[, which(!apply(res, 2, function(x){ all(x == 0) }))]
#dim(res)

op <- par(mfrow=c(4,4))

for(rowm in 1:16){
  plot(x = c(1,max(chrs.starts)), y = c(0,nrow(comparison)), t='n', yaxt='n', ylab="", main=rownames(res)[rowm])
  i <- 1
  for(chr in chrs){
    onChr <- which(markerinfo[,"Chr"] == chr)
    mPos <- as.numeric(markerinfo[onChr, "MapInfo"])
    heights <- res[rowm, rownames(markerinfo[onChr,])]
    points((mPos + chrs.starts[chr]), y = heights, type='l', col=1+((i%%2)==0))
    i <- i + 1
  }
}


for(rowm in 1:16){
  plot(x = c(1,max(chrs.starts)), y = c(0,nrow(comparison)), t='n', yaxt='n', ylab="", main=rownames(res)[rowm])
  i <- 1
  for(chr in chrs){
    onChr <- which(markerinfo[,"Chr"] == chr)
    mPos <- as.numeric(markerinfo[onChr, "MapInfo"])
    heights <- res[rowm, rownames(markerinfo[onChr,])]
    points((mPos + chrs.starts[chr]), y = heights, type='l', col=1+((i%%2)==0))
    i <- i + 1
  }
}

#axis(2, at=1:nrow(comparison), comparison[,2], las=2, cex=0.5)




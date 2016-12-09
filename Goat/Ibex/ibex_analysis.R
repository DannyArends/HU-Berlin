# Analysis of the Ibes Goat SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2016
# first written Nov, 2016

# Load data from Siham
setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
chir1 <- read.csv("FilteredLocationCHIR1.0.txt", sep="\t", row.names=1)
chir1 <- cbind(chir1, Pos = (chir1[,"Start"] + chir1[,"Stop"])/2)
map <- chir1[,c("chrN", "Pos")]

sihamrawsnp  <- read.table("filtered_snps.txt", sep="\t", check.names=FALSE, colClasses="character")
sihamsamples <- read.csv("merged_samples_NO_DN2.txt",sep="\t")

tagData <- rownames(sihamsamples)[which(sihamsamples[,"Breed"] == "Tagg")]
sihamsamples <- sihamsamples[tagData,]
sihamrawsnp <- sihamrawsnp[,tagData]

# Load the Ibex data
setwd("D:/Edrive/Goat/DNA/Ibex")
snpdata <- read.csv("Goat_Nov2016_FinalReport.txt", sep='\t', skip=9, header=TRUE, na.strings=c("","--", "na", "NA"), row.names=1, colClasses="character")
colnames(snpdata) <- paste0("Ind", 1:48)

samples <- read.csv("Goat_Nov2016_Samples.txt", sep="\t", row.names = 1, colClasses="character")
rownames(samples) <- paste0("Ind", 1:48)

snpinfo <- read.csv("Goat_Nov2016_SNP.txt",sep="\t", row.names=1, colClasses="character")

notOnMap <- rownames(snpdata)[which(!(rownames(snpdata) %in% rownames(map)))]
snpdata <- snpdata[-which(rownames(snpdata) %in% notOnMap), ]
snpinfo <- snpinfo[-which(rownames(snpinfo) %in% notOnMap), ]

notOnArray <- rownames(map)[which(!rownames(map) %in% rownames(snpdata))]

map <- map[with(map, order(chrN, Pos)), ]  # Order the map, chromosome then position
snpdata <- snpdata[rownames(map), ]        # Order the SNP data, as to match the map
snpinfo <- snpinfo[rownames(map), ]        # Order the SNP info to match the map

highQualityS <- rownames(samples)[which(as.numeric(samples[, "Call.Rate"]) > .85)]
lowQualityS <- rownames(samples)[which(as.numeric(samples[, "Call.Rate"]) < .60)]

snpdataHQ <- snpdata[, highQualityS]
snpdataLQ <- snpdata[, lowQualityS]

dim(map); map[1:5,]
dim(snpdata); snpdata[1:5, 1:5]
dim(snpinfo); snpinfo[1:5,]

missingHQ <- apply(apply(snpdataHQ, 1, is.na), 2, sum) / ncol(snpdataHQ)

tooMuchMissingHQ <- rownames(snpdataHQ)[which(missingHQ > 0.05)]
tooLowGenTrain   <- rownames(snpinfo)[which(snpinfo[,"GenTrain.Score"] < 0.6)]
tooLowMinorA     <- rownames(snpinfo)[which(snpinfo[,"Minor.Freq"] < 0.05)]

length(tooLowGenTrain)
length(tooMuchMissingHQ) - length(which(tooMuchMissingHQ %in% tooLowGenTrain))
length(tooLowMinorA) - length(which(tooLowMinorA %in% c(tooMuchMissingHQ, tooLowGenTrain)))


badMarkers <- unique(c(tooMuchMissingHQ, tooLowGenTrain, tooLowMinorA))
length(badMarkers)

map <-  map[-which(rownames(map) %in% badMarkers), ]
snpdata <- snpdata[-which(rownames(snpdata) %in% badMarkers), ]
snpdataHQ <- snpdataHQ[-which(rownames(snpdataHQ) %in% badMarkers), ]
snpdataLQ <- snpdataLQ[-which(rownames(snpdataLQ) %in% badMarkers), ]
snpinfo <- snpinfo[-which(rownames(snpinfo) %in% badMarkers), ]

#apply(apply(snpdataLQ, 1, is.na), 2, sum) / ncol(snpdataLQ)

#goodSNPLQ <- rownames(snpdataLQ)[which(apply(apply(snpdataLQ, 1, is.na), 2, sum) / ncol(snpdataLQ) < 0.1)]
#snpdata <- snpdata[goodSNPLQ, ]

#snpdata <- snpdataHQ

dim(map); map[1:5,]
dim(snpdataHQ); snpdataHQ[1:5,1:5]
dim(snpinfo); snpinfo[1:5,]

## Merge the SNP data from Sudanese Ibes with the taggar goat
sihamrawsnp <- sihamrawsnp[which(rownames(sihamrawsnp) %in%  rownames(snpdata)),]
snpdata <- snpdata[which(rownames(snpdata) %in%  rownames(sihamrawsnp)),]

map <- map[rownames(snpdata), ]
map <- map[with(map, order(chrN, Pos)), ]         # Order the map, chromosome then position
snpdata <- snpdata[rownames(map), ]               # Order the SNP data, as to match the map
sihamrawsnp <- sihamrawsnp[rownames(map), ]       # Order the SNP data, as to match the map
snpinfo <- snpinfo[rownames(map), ]               # Order the SNP info to match the map
snpdata <- cbind(snpdata, sihamrawsnp)

for(x in c(1:29, "X")){
  diffs <- diff(map[which(map[,"chrN"] == x), "Pos"])
  cat(x, round(mean(diffs) / 1000,1), round(min(diffs) / 1000,1), "/", round(max(diffs) / 1000,1),"\n")
}

snpinfo <- cbind(snpinfo, minorallele = NA, majorallele = NA)
for(x in rownames(snpdata)){
  rr <- table(unlist(strsplit(unlist(snpdata[x,]), "")))
  snpinfo[x, "minorallele"] <- names(which.min(rr))
  snpinfo[x, "majorallele"] <- names(which.max(rr))
}

numsnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames=list(rownames(snpdata), colnames(snpdata)))
for(x in rownames(snpdata)){
  minor.Homo <- paste0(snpinfo[x, "minorallele"], snpinfo[x, "minorallele"])
  major.Homo <- paste0(snpinfo[x, "majorallele"], snpinfo[x, "majorallele"])
  numsnpdata[x, snpdata[x,] == minor.Homo] <- 1
  numsnpdata[x, snpdata[x,] == major.Homo] <- 3
  numsnpdata[x, snpdata[x,] != major.Homo & snpdata[x,] != minor.Homo & !is.na(snpdata[x,])] <- 2
}
numsnpdataALL <- numsnpdata

numsnpdataHQ <- numsnpdata[,highQualityS]
numsnpdataHQ <- cbind(numsnpdataHQ, numsnpdata[,tagData])

colnames(numsnpdataHQ)[which(grepl("Ind", colnames(numsnpdataHQ)))] <- samples[colnames(numsnpdataHQ)[which(grepl("Ind", colnames(numsnpdataHQ)))], "Sum"]
colnames(numsnpdataHQ)[which(colnames(numsnpdataHQ) %in% tagData)] <- rep("Taggar", length(tagData))


distances <- dist(t(numsnpdataHQ))
plot(hclust(distances), main="Euclidean distances between populations", hang=-1)

## Principal component analysis
misdata <- which(apply(numsnpdataHQ, 1, function(x){any(is.na(x))}))
numsnppca <- numsnpdataHQ[-misdata,]

nonseg <- which(apply(numsnppca, 1, function(x){ length(table(x)) }) == 1 )
numsnppca <- numsnppca[-nonseg,]
dim(numsnppca)

pcares <- prcomp(t(numsnppca), scale=TRUE)
groups <- colnames(numsnppca)

sumpca <- summary(pcares)

pca1 <- paste0("PC1 (", round(sumpca$importance[2,1] * 100,1), "%", " var explained)")
pca2 <- paste0("PC2 (", round(sumpca$importance[2,2] * 100,1), "%", " var explained)")
pca3 <- paste0("PC3 (", round(sumpca$importance[2,3] * 100,1), "%", " var explained)")

# Create colors
types <- c("x","o","#", "t")
names(types) <- c("Ibex (Sudan)", "Ibex (Zoo)", "Bezoar", "Taggar")
# Create colors
cols <- c("red", "blue", "orange", "black")
names(cols) <- c("Ibex (Sudan)", "Ibex (Zoo)", "Bezoar", "Taggar")

#png("PCAplot.png", width=600, height=600, res=300, pointsize = 5)
plot(c(-250,100), c(-100,150), col = cols[as.character(as.character(groups))],pch = 19, xlab=paste0("PCA 1 ",pca1), ylab=paste0("PCA 2 ",pca2), 
      t = 'n',xaxt='n', yaxt='n', main="PCA analysis", cex.lab=0.8)
axis(1, at = seq(-50, 100, 20),cex.axis=0.8)
#abline(v = seq(-50, 100, 15), col="gray", lty=2)
axis(2, at = seq(-100, 150, 20),las=2,cex.axis=0.8)
#abline(h = seq(-100, 150, 50), col="gray", lty=2)
op <- par(mfrow = c(2,2))
plot(pcares$x[,1], pcares$x[,2], col = cols[as.character(as.character(groups))], pch = types[as.character(groups)], cex=1.2, xlab = pca1, ylab = pca2)
plot(pcares$x[,1], pcares$x[,3], col = cols[as.character(as.character(groups))], pch = types[as.character(groups)], cex=1.2, xlab = pca1, ylab = pca3)
plot(pcares$x[,2], pcares$x[,3], col = cols[as.character(as.character(groups))], pch = types[as.character(groups)], cex=1.2, xlab = pca2, ylab = pca3)
plot(c(0,1),c(0,1), t = 'n',xaxt='n',yaxt='n', xlab="", ylab="")
legend("center", unique(groups), pch = types[unique(as.character(groups))], col= cols[unique(as.character(groups))], bg="white")
#dev.off()

 # Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){ return(var.loadings * comp.sdev) }

# Variable correlation/coordinates
var.coord <- t(apply(pcares$rotation, 1, var_cor_func, pcares$sdev))
head(var.coord[, 1:4])
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
colz <- c(rep("black", nrow(var.contrib)-100), rep("red", 100))
op <- par(mfrow = c(2,2))
plot(sort(var.contrib[,1]), main="Sorted contribution to PC1", las=2, xlab="SNP", ylab="% Contribution", col=colz,pch=18)
abline(h=mean(var.contrib[,1]))
abline(h = mean(var.contrib[,1]) + sd(var.contrib[,1]), lty=2)
abline(h = mean(var.contrib[,1]) + 2 * sd(var.contrib[,1]), lty=2, col="green")
abline(h = mean(var.contrib[,1]) - sd(var.contrib[,1]), lty=2)
plot(sort(var.contrib[,2]), main="Sorted contribution to PC2", las=2, xlab="SNP", ylab="% Contribution", col=colz,pch=18)
abline(h=mean(var.contrib[,2]))
abline(h = mean(var.contrib[,2]) + sd(var.contrib[,2]), lty=2)
abline(h = mean(var.contrib[,2]) + 2 * sd(var.contrib[,2]), lty=2, col="green")
abline(h = mean(var.contrib[,2]) - sd(var.contrib[,2]), lty=2)
plot(sort(var.contrib[,3]), main="Sorted contribution to PC3", las=2, xlab="SNP", ylab="% Contribution", col=colz,pch=18)
abline(h=mean(var.contrib[,3]))
abline(h = mean(var.contrib[,3]) + sd(var.contrib[,3]), lty=2)
abline(h = mean(var.contrib[,3]) + 2 * sd(var.contrib[,3]), lty=2, col="green")
abline(h = mean(var.contrib[,3]) - sd(var.contrib[,3]), lty=2)
plot(c(0,1),c(0,1), t = 'n',xaxt='n',yaxt='n', xlab="", ylab="", bty="n", box=FALSE)
legend("center", c("Mean contribution", "Mean +/- sd", "Mean + 2*sd", "Top 100"), lty=c(1, 2, 2,0), pch=c(NA,NA,NA, 18), col=c(1,1,"green", "red"), bg="white")

pc1thr <- mean(var.contrib[,1]) + 2 * sd(var.contrib[,1])
pc2thr <- mean(var.contrib[,2]) + 2 * sd(var.contrib[,2])
pc3thr <- mean(var.contrib[,3]) + 2 * sd(var.contrib[,3])

importantSNPsP1 <-  names(which(var.contrib[,1] >= pc1thr)); length(importantSNPsP1)
importantSNPsP1t <-  names(sort(var.contrib[,1],dec=TRUE))[1:100]; length(importantSNPsP1t)
importantSNPsP2 <-  names(which(var.contrib[,2] >= pc2thr)); length(importantSNPsP2)
importantSNPsP2t <-  names(sort(var.contrib[,2],dec=TRUE))[1:100]; length(importantSNPsP2t)
importantSNPsP3 <-  names(which(var.contrib[,3] >= pc3thr)); length(importantSNPsP3)
importantSNPsP3t <-  names(sort(var.contrib[,3],dec=TRUE))[1:100]; length(importantSNPsP3t)

SNPsPCA1 <- map[importantSNPsP1, c("chrN", "Pos")]
SNPsPCA1t <- map[importantSNPsP1t, c("chrN", "Pos")]
SNPsPCA2 <- map[importantSNPsP2, c("chrN", "Pos")]
SNPsPCA2t <- map[importantSNPsP2t, c("chrN", "Pos")]
SNPsPCA3 <- map[importantSNPsP3, c("chrN", "Pos")]
SNPsPCA3t <- map[importantSNPsP3t, c("chrN", "Pos")]

chromosomes <- c(as.character(1:29),"X")
ymax <- max(map[,"Pos"])

op <- par(mfrow = c(1,1))
plot(x=c(0,length(chromosomes)), y = c(0, ymax), t='n', xaxt='n', yaxt='n', ylab="", xlab="Chromosome")
chrid <- 1
for(chr in chromosomes){
  allG <- map[map[,"chrN"] == chr, "Pos"]
  chrM <- max(map[map[,"chrN"] == chr, "Pos"])
  intG1 <- SNPsPCA1[which(SNPsPCA1[,"chrN"] == chr), "Pos"]
  intG1t <- SNPsPCA1t[which(SNPsPCA1t[,"chrN"] == chr), "Pos"]
  intG2 <- SNPsPCA2[which(SNPsPCA2[,"chrN"] == chr), "Pos"]
  intG2t <- SNPsPCA2t[which(SNPsPCA2t[,"chrN"] == chr), "Pos"]
  intG3 <- SNPsPCA3[which(SNPsPCA3[,"chrN"] == chr), "Pos"]
  intG3t <- SNPsPCA3t[which(SNPsPCA3t[,"chrN"] == chr), "Pos"]
  lines(x=c(chrid,chrid), y = c(0, chrM),lwd=2)
  #points(x = rep(chrid, length(allG))+0.0, allG, pch="-", col="lightgray", cex=0.5)
  points(x = rep(chrid, length(intG1))-0.2, intG1, pch="-", col=rgb(1,0,0,0.1), cex=2)
  points(x = rep(chrid, length(intG1t))-0.2, intG1t, pch="-", col="red", cex=2)
  points(x = rep(chrid, length(intG2))+0.0, intG2, pch="-", col="gray90", cex=2)
  points(x = rep(chrid, length(intG2t))+0.0, intG2t, pch="-", col="black", cex=2)
  points(x = rep(chrid, length(intG3))+0.2, intG3, pch="-", col="aliceblue", cex=2)
  points(x = rep(chrid, length(intG3t))+0.2, intG3t, pch="-", col="blue", cex=2)
  chrid <- chrid + 1
}
axis(1, at = 1:length(chromosomes), chromosomes,cex.axis=0.8)
axis(2, at = seq(0, ymax, 10000000), paste(seq(0, ymax, 10000000) / 1000000, "mb"),las=2)
legend("top", c("PC1 > (mean + 2SD)", "PC1 - Top 100", "PC2 > (mean + 2SD)", "PC2 - Top 100", "PC3 > (mean + 2SD)", "PC3 - Top 100"), pch="-", col=c(rgb(1,0,0,0.1), "red", "gray90", "black", "aliceblue", "blue"),pt.cex=3)

pc1map <- map[rownames(SNPsPCA1t),]
pc2map <- map[rownames(SNPsPCA2t),]
pc3map <- map[rownames(SNPsPCA3t),]

proteins <- read.csv("D:/Edrive/Goat/DNA/annotation/ProteinTable10731_39633_ensembl.txt", sep="\t")

getGenes <- function(pcmap, proteins){
  genes <- NULL
  for(x in 1:nrow(pcmap)){
    inregion <- proteins[which(as.character(proteins[,1]) == pcmap[x,"chrN"] & 
                               proteins[,3] > as.numeric(pcmap[x,"Pos"]) - 250000 & 
                               proteins[,3] < as.numeric(pcmap[x,"Pos"]) + 250000),]
    genes <- rbind(genes, inregion)
  }
  gn <- unlist(lapply(strsplit(as.character(genes[,"Protein.name"]), " isoform"),"[",1))
  genes <- genes[-which(duplicated(gn)),]
  genes <- genes[with(genes, order(X.Replicon.Name, Start)), ]         # Order the map, chromosome then position
  return(genes)
}

write.table(getGenes(pc1map, proteins), file="pc1genes.txt", sep="\t")
write.table(getGenes(pc2map, proteins), file="pc2genes.txt", sep="\t")
write.table(getGenes(pc3map, proteins), file="pc3genes.txt", sep="\t")

## OLD CLUSTERING PLOTS
distances <- dist(t(numsnpdata))
plot(hclust(distances), main="Clustering of High Quality data")

distances <- dist(t(numsnpdata[, highQualityS]))
plot(hclust(distances), main="Clustering of High Quality data")


snpdata <- snpdata[, c(highQualityS, tagData)]

# Transform SNPs to a format stampp can understand
absnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames=list(rownames(snpdata), colnames(snpdata)))
for(x in rownames(snpdata)){
  minor.Homo <- paste0(snpinfo[x, "minorallele"], snpinfo[x, "minorallele"])
  major.Homo <- paste0(snpinfo[x, "majorallele"], snpinfo[x, "majorallele"])
  absnpdata[x, snpdata[x,] == minor.Homo] <- "AA"
  absnpdata[x, snpdata[x,] == major.Homo] <- "BB"
  absnpdata[x, snpdata[x,] != major.Homo & snpdata[x,] != minor.Homo & !is.na(snpdata[x,])] <- "AB"
}
write.table(absnpdata, "absnpdata.txt", sep="\t", quote=FALSE)

library(StAMPP)

stammpinput <- t(absnpdata)
populations <- as.character(samples[rownames(stammpinput), "Sum"])
populations[is.na(populations)] <- "Taggar"

stammpinput <- cbind(Sample = rownames(stammpinput), Pop = populations, Ploidy = 2, Format = "BiA", stammpinput)
stammpinput <- as.data.frame(stammpinput)

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values
stammp.D.ind <- stamppNeisD(stammpinput.freq, FALSE) # Population D values
stammpinput.fst <- stamppFst(stammpinput.freq, 1000, 95, 4) # Population Fst values
stammpinput.fst$Fsts
write.table(stammpinput.fst$Fsts, file = "fstsHQ.txt", sep = "\t")

stammpinput.amova <- stamppAmova(stammp.D.ind, stammpinput.freq, 10000)
write.table(stammpinput.amova[[1]], file = "amovaSSD.txt", sep = "\t")
write.table(stammpinput.amova[[3]], file = "amovapvalues.txt", sep = "\t")

# Nei's genetic distance
rownames(stammp.D.ind) <- samples[rownames(stammp.D.ind), "Sum"]
rownames(stammp.D.ind)[is.na(rownames(stammp.D.ind))] <- "Taggar"
rownames(stammp.D.ind)[is.na(rownames(stammp.D.ind))] <- "Taggar"
stmpD <- as.dist(stammp.D.ind)
plot(hclust(stmpD), main="Nei's genetic distance")


### diversity analysis
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
genotypes_genpopA <- t(toGenPop(absnpdata))

write.table(genotypes_genpopA, "snpdataGP.txt", sep="\t", quote=FALSE)
set.seed(1)
rsample <- sample(ncol(genotypes_genpopA), 2000)
genotypes_genpop <- genotypes_genpopA[, rsample]

rownames(genotypes_genpop) <- samples[rownames(genotypes_genpop), "Sum"]
rownames(genotypes_genpop)[is.na(rownames(genotypes_genpop))] <- "Taggar"
rownames(genotypes_genpop) <- gsub(")", "", gsub("(", "", gsub(" ", "_", rownames(genotypes_genpop)), fixed = TRUE))

### Write out genotypes in genpop format for usage in diveRsity
cat("BLANK\n", file="genotypes_genpop.txt")
cat(paste0(colnames(genotypes_genpop), collapse="\n"), file="genotypes_genpop.txt", append=TRUE)
cat("\n", file="genotypes_genpop.txt", append=TRUE)

for(pop in unique(rownames(genotypes_genpop))) {
  cat("POP\n", file="genotypes_genpop.txt", append=TRUE)
  ii <- which(rownames(genotypes_genpop) == pop)
  write.table(genotypes_genpop[ii,], sep = " ", quote=FALSE, na = "0000", file="genotypes_genpop.txt", append=TRUE, col.names=FALSE)
}
cat("POP\n", file="genotypes_genpop.txt", append=TRUE)
genotypes_copy <- genotypes_genpop
rownames(genotypes_copy) <- rep("Combined", nrow(genotypes_copy))
write.table(genotypes_copy[,], sep = " ", quote=FALSE, na = "0000", file="genotypes_genpop.txt", append=TRUE, col.names=FALSE)

library(diveRsity)
basicStats <- divBasic(infile = "genotypes_genpop.txt", outfile="basicStatsOut.txt", gp=2, bootstraps = 100)
names(basicStats$fis) <- colnames(basicStats$Ho)

basicStats$Ho["overall",] # Observed
basicStats$He["overall",] # Expected

basicStats$fis[["Bezoar"]]["overall",]                   # Fis
basicStats$fis[["Ibex_Zoo"]]["overall",]                      # Fis
basicStats$fis[["Ibex_Sudan"]]["overall",]                         # Fis
basicStats$fis[["Combined"]]["overall",]                      # Fis

advancedStats <- diffCalc(infile = "genotypes_genpop.txt", outfile="fstOnlyOut.txt", fst=TRUE, pairwise=TRUE, boots = 1000)
advancedStats$pairwise$Fst # Pairwise Fst per populations

## STRUCTURE
if(!file.exists("cleaned_genotypes_structure.txt")){
  # Write out the data for STRUCTURE
  numsnpdata <- t(numsnpdataALL[, c(highQualityS, tagData)])
  structGeno <- NULL #matrix(NA, nrow(numGeno) * 2, ncol(numGeno))
  for(x in 1:nrow(numsnpdata)){
    gg <- rbind(rep(NA, ncol(numsnpdata)), rep(NA, ncol(numsnpdata)))
    a1 <- which(numsnpdata[x,] == 1)
    a2 <- which(numsnpdata[x,] == 2)
    a3 <- which(numsnpdata[x,] == 3)
    gg[1, a1] <- 0; gg[2, a1] <- 0  # Pretty inefficient, but it will do the trick
    gg[1, a2] <- 0; gg[2, a2] <- 1
    gg[1, a3] <- 1; gg[2, a3] <- 1
    gg[is.na(gg)] <- 9
    structGeno <- rbind(structGeno, gg)
  }
 
  rownames(structGeno) <- gsub(" ","", unlist(lapply(rownames(numsnpdata), rep, 2)))    # Spaces are not allowed in sample names
  colnames(structGeno) <- colnames(numsnpdata)

  rownames(samples) <- gsub(" ","", rownames(samples))
  samples[gsub(" ","",rownames(structGeno)),]
  
  pops <- samples[gsub(" ","",rownames(structGeno)),"Sum"]
  pops[is.na(pops)] <- "Taggar"
  pops <- as.factor(pops)
  structGeno <- cbind(as.numeric(pops), structGeno)
  # We need to remove the empty entry from the first line of the file
  # Structure settings: 44 individuals, 33698 markers, missing data = 9, row markernames: yes, column individual id: yes, column pop id: yes
  write.table(structGeno, file="cleaned_genotypes_structure.txt", sep = "\t")           # Save the genotypes to disk
}


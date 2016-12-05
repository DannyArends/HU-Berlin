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

snpdata <- cbind(snpdata, sihamrawsnp[,tagData])

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
plot(hclust(distances), main="Clustering of High Quality data")

## Principal component analysis
misdata <- which(apply(numsnpdataHQ, 1, function(x){any(is.na(x))}))
numsnppca <- numsnpdataHQ[-misdata,]

nonseg <- which(apply(numsnppca, 1, function(x){ length(table(x)) }) == 1 )
numsnppca <- numsnppca[-nonseg,]
dim(numsnppca)

pcares <- prcomp(t(numsnppca), scale=TRUE)
groups <- colnames(numsnppca)

sumpca <- summary(pcares)

pca1 <- paste0("(", round(sumpca$importance[2,1] * 100,1), "%", " var explained)")
pca2 <- paste0("(", round(sumpca$importance[2,2] * 100,1), "%", " var explained)")

# Create colors
types <- c("x","o","#", "t")
names(types) <- c("Sudan", "Zoo Ibex", "Bezoarziege", "Taggar")
# Create colors
cols <- c("red", "blue", "orange", "black")
names(cols) <- c("Sudan", "Zoo Ibex", "Bezoarziege", "Taggar")

#png("PCAplot.png", width=600, height=600, res=300, pointsize = 5)
plot(c(-250,100), c(-100,150), col = cols[as.character(as.character(groups))],pch = 19, xlab=paste0("PCA 1 ",pca1), ylab=paste0("PCA 2 ",pca2), 
      t = 'n',xaxt='n', yaxt='n', main="PCA analysis", cex.lab=0.8)
axis(1, at = seq(-50, 100, 20),cex.axis=0.8)
#abline(v = seq(-50, 100, 15), col="gray", lty=2)
axis(2, at = seq(-100, 150, 20),las=2,cex.axis=0.8)
#abline(h = seq(-100, 150, 50), col="gray", lty=2)
plot(pcares$x[,1], pcares$x[,2], col = cols[as.character(as.character(groups))], pch = types[as.character(groups)], cex=0.6)
#legend("topright", c("Taggar", "Desert", "Nilotic", "Nubian"), col=cols, pch=types, bg="white", cex=0.8)
#dev.off()

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
importantSNPsP1 <-  names(which(var.contrib[,1] >= 0.0083))
#importantSNPsP2sh <-  names(which(var.contrib[,2] >= 0.05))
#importantSNPsP2h <-  names(which(var.contrib[,2] >= 0.04))
SNPsPCA1 <- map[importantSNPsP1, c("chrN", "Pos")]
#SNPsPCA2 <- map[importantSNPsP2h, c("chrN", "Pos")]

chromosomes <- c(as.character(1:29),"X")
ymax <- max(map[,"Pos"])

plot(x=c(0,length(chromosomes)), y = c(0, ymax), t='n', xaxt='n', yaxt='n', ylab="", xlab="Chromosome")
chrid <- 1
for(chr in chromosomes){
  allG <- map[map[,"chrN"] == chr, "Pos"]
  chrM <- max(map[map[,"chrN"] == chr, "Pos"])
  intG1 <- SNPsPCA1[which(SNPsPCA1[,"chrN"] == chr),"Pos"]
  #intG2 <- SNPsPCA2[which(SNPsPCA2[,"chrN"] == chr),"Pos"]
  #intP <- SNPlPCA1[which(SNPlPCA1[,"chrN"] == chr),"Pos"]
  lines(x=c(chrid,chrid), y = c(0, chrM))
  points(x = rep(chrid, length(allG)), allG, pch="-", col=rgb(0.9, 0.9, 0.9), cex=2)
  #points(x = rep(chrid, length(intP)), intP, pch="-", col=rgb(0.5, 0.5, 0.5), cex=2)
  #points(x = rep(chrid, length(intG2)), intG2, pch="-", col="darkgray", cex=2)
  points(x = rep(chrid, length(intG1)), intG1, pch="-", col="orange", cex=2)
  chrid <- chrid + 1
}
axis(1, at = 1:length(chromosomes), chromosomes,cex.axis=0.8)
axis(2, at = seq(0, ymax, 10000000), paste(seq(0, ymax, 10000000) / 1000000, "mb"),las=2)

HighIbexAll <- numsnpdataALL[rownames(SNPsPCA1),]
colnames(HighIbexAll) <- samples[colnames(numsnpdataALL[rownames(SNPsPCA1),]),"Sum"]

highBetweenIbex <- numsnppca[rownames(SNPsPCA1),]
#colnames(highBetweenIbex) <- samples[colnames(numsnppca[rownames(SNPsPCA1),]),"Sum"]

res <- map[rownames(highBetweenIbex),]

proteins <- read.csv("D:/Edrive/Goat/DNA/annotation/ProteinTable10731_39633_ensembl.txt", sep="\t")

genes <- NULL
for(x in 1:nrow(res)){
  inregion <- proteins[which(as.character(proteins[,1]) == res[x,"chrN"] & 
                             proteins[,3] > as.numeric(res[x,"Pos"]) - 250000 & 
                             proteins[,3] < as.numeric(res[x,"Pos"]) + 250000),]
  genes <- rbind(genes, inregion)
}
genes <- genes[-which(duplicated(genes[,"GeneID"])),]

write.table(genes, file="res.txt", sep="\t")


#colnames(numsnpdata) <- paste0(samples[colnames(snpdata),"Origin.Species"]," ", round(100 * (1-as.numeric(samples[colnames(snpdata),"Call.Rate"])), d = 0))
#idx <- grep("Sudan", colnames(numsnpdata))
#numsnpdata <- numsnpdata[, idx]

distances <- dist(t(numsnpdata))
plot(hclust(distances), main="Clustering of High Quality data")


distances <- dist(t(numsnpdata[, highQualityS]))
plot(hclust(distances), main="Clustering of High Quality data")


# Transform SNPs to a format stampp can understand
absnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames=list(rownames(snpdata), colnames(snpdata)))
for(x in rownames(snpdata)){
  minor.Homo <- paste0(snpinfo[x, "minorallele"], snpinfo[x, "minorallele"])
  major.Homo <- paste0(snpinfo[x, "majorallele"], snpinfo[x, "majorallele"])
  absnpdata[x, snpdata[x,] == minor.Homo] <- "AA"
  absnpdata[x, snpdata[x,] == major.Homo] <- "BB"
  absnpdata[x, snpdata[x,] != major.Homo & snpdata[x,] != minor.Homo & !is.na(snpdata[x,])] <- "AB"
}

library(StAMPP)

stammpinput <- t(absnpdata)
stammpinput <- cbind(Sample = rownames(stammpinput), Pop = as.character(samples[rownames(stammpinput), "Sum"]), Ploidy = 2, Format = "BiA", stammpinput)
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
rownames(stammp.D.ind) <- samples[rownames(stammp.D.ind), "Origin.Species"]
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

genotypes_genpop <- t(toGenPop(absnpdata))
set.seed(0)
rsample <- sample(ncol(genotypes_genpop), ncol(genotypes_genpop) / 20)
genotypes_genpop <- genotypes_genpop[, rsample]

rownames(genotypes_genpop) <- samples[rownames(genotypes_genpop), "Sum"]



### Write out genotypes in genpop format for usage in diveRsity
cat("BLANK\n", file="genotypes_genpop.txt")
cat(paste0(colnames(genotypes_genpop), collapse="\n"), file="genotypes_genpop.txt", append=TRUE)

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

basicStats$fis[["Bezoarziege"]]["overall",]                   # Fis
basicStats$fis[["Zoo Ibex"]]["overall",]                      # Fis
basicStats$fis[["Sudan"]]["overall",]                         # Fis
basicStats$fis[["Combined"]]["overall",]                      # Fis


advancedStats <- diffCalc(infile = "genotypes_genpop.txt", outfile="fstOnlyOut.txt", fst=TRUE, pairwise=TRUE, boots = 1000)
advancedStats$pairwise$Fst # Pairwise Fst per populations

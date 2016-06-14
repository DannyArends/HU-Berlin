# Analysis of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016
library(ape)

setwd("E:/Goat/DNA/SihamAnalysis")

snpdata    <- read.table("filtered_snps.txt", sep="\t", check.names=FALSE, colClasses="character")
snpinfo    <- read.table("snpinfo.txt", sep="\t", na.strings=c("", "NA", "N.D."))
samples    <- read.table("sampleinfo.txt", sep="\t")
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`
samples    <- cbind(samples, locations[rownames(samples),])
samples    <- cbind(samples, locationShort = as.character(unlist(lapply(strsplit(as.character(samples[,"Location"]), "_"), "[",1))))

snpdata <- snpdata[,-which(colnames(snpdata) == "DN 2")] # Throw away the duplicate individual because it confuses STRUCTURE
samples <- samples[-which(rownames(samples) == "DN 2"),] # Throw away the duplicate individual because it confuses STRUCTURE

write.table(samples, "merged_samples_NO_DN2.txt",sep="\t")

snpAlleles <- lapply(strsplit(as.character(snpinfo[,"allele"]), ""), "[", c(1,3))

chir1 <- read.csv("FilteredLocationCHIR1.0.txt", sep="\t", row.names=1)
chir1 <- cbind(chir1, Pos = (chir1[,"Start"] + chir1[,"Stop"])/2)
chir2 <- read.csv("FilteredLocationCHIR2.0.txt", sep="\t", row.names=1)
chir2 <- cbind(chir2, Pos = (chir2[,"Start"] + chir2[,"Stop"])/2)

snpinfo <- cbind(snpinfo, Chr_C1 = NA)
snpinfo <- cbind(snpinfo, Pos_C1 = NA)

snpinfo <- cbind(snpinfo, Chr_C2 = NA)
snpinfo <- cbind(snpinfo, Pos_C2 = NA)

snpinfo[rownames(snpinfo), "Chr_C1"] <- chir1[rownames(snpinfo),"chrN"]
snpinfo[rownames(snpinfo), "Pos_C1"] <- chir1[rownames(snpinfo),"Pos"]

snpinfo[rownames(snpinfo), "Chr_C2"] <- chir2[rownames(snpinfo),"chrN"]
snpinfo[rownames(snpinfo), "Pos_C2"] <- chir2[rownames(snpinfo),"Pos"]

write.table(snpinfo, "merged_snp_info.txt",sep="\t")

if(!file.exists("filtered_snps_numeric_NO_DN2.txt")){
  numsnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames = list(rownames(snpdata), colnames(snpdata)))
  for(x in 1:length(snpAlleles)) {
    if(!is.na(snpinfo[x, "reference"]) && snpAlleles[[x]][1] !=  snpinfo[x, "reference"]){  # C/T while reference is T, so flip it around
      snpAlleles[[x]] <- snpAlleles[[x]][2:1]
    }

    g1 <- paste(snpAlleles[[x]][1], snpAlleles[[x]][1],sep="")
    g2a <- paste(snpAlleles[[x]][1], snpAlleles[[x]][2],sep="")
    g2b <- paste(snpAlleles[[x]][2], snpAlleles[[x]][1],sep="")
    g3 <- paste(snpAlleles[[x]][2], snpAlleles[[x]][2],sep="")
    if(!all(snpdata[x,] %in% c(g1,g2a,g2b,g3, NA))) stop("Nope")
    numsnpdata[x, which(snpdata[x, ] == g1)] <- 1
    numsnpdata[x, which(snpdata[x, ] == g2a)] <- 2
    numsnpdata[x, which(snpdata[x, ] == g2b)] <- 2
    numsnpdata[x, which(snpdata[x, ] == g3)] <- 3
  }

  write.table(numsnpdata, "filtered_snps_numeric_NO_DN2.txt", sep="\t", quote=FALSE)
}else{
  numsnpdata <- read.csv("filtered_snps_numeric_NO_DN2.txt", sep="\t", check.names=FALSE)
}

if(!file.exists("filtered_snps_AB_NO_DN2.txt")){
  absnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames = list(rownames(snpdata), colnames(snpdata)))
  for(x in 1:length(snpAlleles)) {
    if(!is.na(snpinfo[x, "reference"]) && snpAlleles[[x]][1] !=  snpinfo[x, "reference"]){  # C/T while reference is T, so flip it around
      snpAlleles[[x]] <- snpAlleles[[x]][2:1]
    }

    g1 <- paste(snpAlleles[[x]][1], snpAlleles[[x]][1],sep="")
    g2a <- paste(snpAlleles[[x]][1], snpAlleles[[x]][2],sep="")
    g2b <- paste(snpAlleles[[x]][2], snpAlleles[[x]][1],sep="")
    g3 <- paste(snpAlleles[[x]][2], snpAlleles[[x]][2],sep="")
    if(!all(snpdata[x,] %in% c(g1,g2a,g2b,g3, NA))) stop("Nope")
    absnpdata[x, which(snpdata[x, ] == g1)]  <- "AA"
    absnpdata[x, which(snpdata[x, ] == g2a)] <- "AB"
    absnpdata[x, which(snpdata[x, ] == g2b)] <- "AB"
    absnpdata[x, which(snpdata[x, ] == g3)]  <- "BB"
  }

  write.table(absnpdata, "filtered_snps_AB_NO_DN2.txt", sep="\t", quote=FALSE)
}else{
  absnpdata <- read.csv("filtered_snps_AB_NO_DN2.txt", sep="\t", check.names=FALSE)
}

breeds <- as.character(unique(samples[,"Breed"]))
# Minor Allele Frequencies for different breeds
MAFs <- matrix(NA, nrow(numsnpdata),length(breeds), dimnames=list(rownames(numsnpdata), breeds))
for(breed in breeds){
  individuals <-  rownames(samples)[which(samples[,"Breed"] == breed)]
  MAFs[, breed] <- apply(numsnpdata[, individuals], 1, function(x){
    tabulated <- table(unlist(x))
    ref <- sum(tabulated["1"] * 2, tabulated["2"],na.rm=TRUE)
    alt <- sum(tabulated["3"] * 2, tabulated["2"],na.rm=TRUE)
    if(is.na(ref) || is.na(alt)) return(0)
    if(ref < alt) return(ref / (ref+alt))
    if(ref >= alt) return(alt / (ref+alt))
  })
}

### (Non-)Polymorphic loci per group
apply(MAFs,2, function(x){return(length(which(x < 0.05)))})
apply(MAFs,2, function(x){return(length(which(x >= 0.05)))})

breedSpecific <- names(which(apply(MAFs,1,function(x){(length(which(x == 0)) == 3)})))

## STRUCTURE
if(!file.exists("cleaned_genotypes_structure_NO_DN2.txt")){
  # Write out the data for STRUCTURE
  numsnpdata <- t(numsnpdata)
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
  write.table(structGeno, file="cleaned_genotypes_structure_NO_DN2.txt", sep = "\t")           # Save the genotypes to disk
}

library(StAMPP)

stammpinput <- t(absnpdata)
stammpinput <- cbind(Sample = rownames(stammpinput), Pop = as.character(samples[rownames(stammpinput),"Breed"]), Ploidy = 2, Format = "BiA", stammpinput)
stammpinput <- as.data.frame(stammpinput)

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values
stammp.D.ind <- stamppNeisD(stammpinput.freq, FALSE) # Population D values
stammpinput.fst <- stamppFst(stammpinput.freq, 1000, 95, 4) # Population Fst values
stammpinput.fst$Fsts
write.table(stammpinput.fst$Fsts, file = "fsts.txt", sep = "\t")

stammpinput.amova <- stamppAmova(stammp.D.ind, stammpinput.freq, 100)
write.table(stammpinput.amova[[1]], file = "amovaSSD.txt", sep = "\t")
write.table(stammpinput.amova[[3]], file = "amovapvalues.txt", sep = "\t")

numsnpdata <- t(numsnpdata) # Change it back from the STRUCTURE 

tagg <- rownames(samples)[which(samples[,"Breed"] == "Tagg")]
dese <- rownames(samples)[which(samples[,"Breed"] == "Dese")]
ni <- rownames(samples)[which(samples[,"Breed"] == "Ni")]
nu <- rownames(samples)[which(samples[,"Breed"] == "Nu")]

taggtbl <- apply(numsnpdata[,tagg], 1,table)
taggval <- unlist(lapply(taggtbl, function(x){ return(max(x) / sum(x)) }))
plot(taggval)


### Clustering of data

numsnpclustering <- numsnpdata
colnames(numsnpclustering) <- samples[colnames(numsnpdata), "Breed"]

rownames(stammp.D.ind) <- samples[rownames(stammp.D.ind), "Breed"]
colnames(stammp.D.ind) <- samples[colnames(stammp.D.ind), "Breed"]

# Add the reference animal to the dataset, since reference is always coded as 1
#numsnpdata <- cbind(numsnpdata, Reference = 1)
differences <- dist(t(numsnpclustering), method = "manhattan")
stmpD <- as.dist(stammp.D.ind)

clustering1 <- hclust(differences)
clustering2 <- hclust(stmpD)
#plot(as.phylo(clustering), type="r")  # Or a rooted dendrogram plot: plot(root(as.phylo(clustering),"Reference"), type="r")

breeds <- samples[,"Breed"]
names(breeds) <- rownames(samples)

# Create colors
cols <- c("red", "blue", "orange", "black")
names(cols) <- c("Tagg", "Dese", "Ni", "Nu")

labelCol <- function(x) {
  if (is.leaf(x)) {
    hclass <- as.character(attr(x, "label"))             # Fetch the class label
    hcol <- cols[hclass]                            # Determine color of the label
    cat(attr(x, "label"), hclass, hcol, "\n")
    attr(x, "nodePar") <- list(lab.col=hcol, cex=0.1)
  }
  return(x)
}
dendrogram1 <- as.dendrogram(clustering1)
dendrogram2 <- as.dendrogram(clustering2)
dendrogram1.col <- dendrapply(dendrogram1, labelCol)
dendrogram2.col <- dendrapply(dendrogram2, labelCol)

op <- par(mfrow=c(1,2), cex=0.5)
plot(dendrogram1.col, main = "Manhattan distance",cex.axis=1.4, cex.main=1.4)
plot(dendrogram2.col, main = "Nei's genetic distance",cex.axis=1.4, cex.main=1.4)

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

rownames(genotypes_genpop) <- gsub(" ","", paste0(breeds[rownames(genotypes_genpop)],","))


### Write out genotypes in genpop format for usage in diveRsity

cat("BLANK\n", file="genotypes_genpop.txt")
for(x in colnames(genotypes_genpop)){
  cat(paste0(x, "\n"), file="genotypes_genpop.txt", append=TRUE)
}
for(pop in unique(rownames(genotypes_genpop))) {
  cat("POP\n", file="genotypes_genpop.txt", append=TRUE)
  ii <- which(rownames(genotypes_genpop) == pop)
  write.table(genotypes_genpop[ii,], sep = " ", quote=FALSE, na = "0000", file="genotypes_genpop.txt", append=TRUE, col.names=FALSE)
}
cat("POP\n", file="genotypes_genpop.txt", append=TRUE)
genotypes_copy <- genotypes_genpop
rownames(genotypes_copy) <- rep("Combined", nrow(genotypes_copy))
write.table(genotypes_copy[,], sep = " ", quote=FALSE, na = "0000", file="genotypes_genpop.txt", append=TRUE, col.names=FALSE)

### Diversity analysis

library(diveRsity)
starttime <-  proc.time()

basicStats <- divBasic(infile = "genotypes_genpop.txt", outfile="fstOnlyOut.txt", gp=2, bootstraps = 100)
names(basicStats$fis) <- colnames(basicStats$Ho)

endtime <-  proc.time()

basicStats$Ho["overall",] # Observed
basicStats$He["overall",] # Expected
basicStats$fis[["Tagg,"]]["overall",]              # Fis
basicStats$fis[["Ni,"]]["overall",]              # Fis
basicStats$fis[["Nu,"]]["overall",]              # Fis
basicStats$fis[["Dese,"]]["overall",]              # Fis
basicStats$fis[["Combined"]]["overall",]     # Fis

advancedStats <- diffCalc(infile = "genotypes_genpop.txt", outfile="fstOnlyOut.txt", fst=TRUE, pairwise=TRUE, boots = 1000)
# Pairwise Fst per populations
advancedStats$pairwise$Fst


### Principal component analysis
misdata <- which(apply(numsnpdata, 1, function(x){any(is.na(x))}))
numsnppca <- numsnpdata[-misdata,]

pcares <- prcomp(t(numsnppca), scale=TRUE)
groups <- samples[colnames(numsnppca), "Breed"]

sumpca <- summary(pcares)

pca1 <- paste0("(", round(sumpca$importance[2,1] * 100,1), "%", " var explained)")
pca2 <- paste0("(", round(sumpca$importance[2,2] * 100,1), "%", " var explained)")

plot(c(-50,100), c(-100,150), col = cols[as.character(groups)],pch = 19, xlab=paste0("PCA 1 ",pca1), ylab=paste0("PCA 2 ",pca2), t = 'n',xaxt='n',yaxt='n', main="PCA analysis")
axis(1, at = seq(-50, 100, 5))
#abline(v = seq(-50, 100, 15), col="gray", lty=2)
axis(2, at = seq(-100, 150, 20),las=2)
#abline(h = seq(-100, 150, 50), col="gray", lty=2)
points(pcares$x[,1], pcares$x[,2], col = cols[as.character(groups)], pch = 20, cex=1.4)
legend("topright", names(cols), col=cols, pch=20,bg="white")


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

library(devtools)
library(ggbiplot)
g <- ggbiplot(pcares, choices = c(1, 2), obs.scale = 1, var.scale = 1, groups = groups, ellipse = TRUE, circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)

### Look into the rotations
importantSNPs <-  names(which(var.contrib[,1] >= 0.02))
lessimpSNPs <-  names(which(var.contrib[,1] >= 0.02))

#snpinfo <- snpinfo[-which(!rownames(snpinfo) %in% rownames(pcares$rotation)), ]

SNPsPCA1 <- snpinfo[importantSNPs, c("Chr", "Position")]
SNPlPCA1 <- snpinfo[lessimpSNPs, c("Chr", "Position")]

chromosomes <- c(as.character(1:29),"X")
ymax <- max(snpinfo[,"Position"])

plot(x=c(0,length(chromosomes)), y = c(0, ymax), t='n', xaxt='n', yaxt='n', ylab="", xlab="Chromosome")
chrid <- 1
for(chr in chromosomes){
  allG <- snpinfo[snpinfo[,"Chr"] == chr, "Position"]
  chrM <- max(snpinfo[snpinfo[,"Chr"] == chr, "Position"])
  intG <- SNPsPCA1[which(SNPsPCA1[,"Chr"] == chr),"Position"]
  intP <- SNPlPCA1[which(SNPlPCA1[,"Chr"] == chr),"Position"]
  lines(x=c(chrid,chrid), y = c(0, chrM))
  points(x = rep(chrid, length(allG)), allG, pch="-", col=rgb(0.9, 0.9, 0.9), cex=2)
  points(x = rep(chrid, length(intP)), intP, pch="-", col=rgb(0.5, 0.5, 0.5), cex=2)
  points(x = rep(chrid, length(intG)), intG, pch="-", col="black", cex=2)
  chrid <- chrid + 1
}
axis(1, at = 1:length(chromosomes), chromosomes,cex.axis=0.8)
axis(2, at = seq(0, ymax, 10000000), paste(seq(0, ymax, 10000000) / 1000000, "mb"),las=2)

SNPsPCA1 <- cbind(SNPsPCA1, var.contrib = var.contrib[row.names(SNPsPCA1),1])

res <- NULL
for(x in chromosomes){
  res <- rbind(res, SNPsPCA1[which(SNPsPCA1[,"Chr"] == x),])
}

proteins <- read.csv("../annotation/ProteinTable10731_39633_ensembl.txt",sep="\t")

genes <- NULL
for(x in 1:nrow(res)){
  inregion <- proteins[which(as.character(proteins[,1]) == res[x,"Chr"] & 
                             proteins[,3] > as.numeric(res[x,"Position"]) - 250000 & 
                             proteins[,3] < as.numeric(res[x,"Position"]) + 250000),]
  genes <- rbind(genes, inregion)
}
genes <- genes[-which(duplicated(genes[,"GeneID"])),]
write.table(genes[,-c(2,7)], file="genesNearbyPCA1snps.txt", sep="\t", row.names=FALSE, quote=FALSE)

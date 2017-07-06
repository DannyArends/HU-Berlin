# Analysis of STRUCTURE results
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jul, 2016
# first written Jul, 2016

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
samples <- read.table("sampleinfo.txt", sep="\t")
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`
samples    <- cbind(samples, locations[rownames(samples),])
samples    <- cbind(samples, locationShort = as.character(unlist(lapply(strsplit(as.character(samples[,"Location"]), "_"), "[",1))))
samples    <- samples[-which(rownames(samples) == "DN 2"),] # Throw away the duplicate individual because it confuses STRUCTURE

popinfo <- as.character(samples[,"Breed"])
names(popinfo) <- gsub(" ", "", rownames(samples))

popinfo[popinfo == "Tagg"] <- "T"
popinfo[popinfo == "Dese"] <- "D"

plotStructure <- function(stmatrix, popinfo, doSort = FALSE, sortTwice = FALSE){
  breaks <- 0
  if(doSort){
    ordering <- NULL
    stmatrix.copy <- stmatrix
    for(x in 1:ncol(stmatrix.copy)){
      sorted <- sort(stmatrix.copy[,x], decreasing = TRUE)
      samples <- names(sorted[sorted >= 0.5])
      ordering <- c(ordering, samples)
      stmatrix.copy <- stmatrix.copy[-which(rownames(stmatrix.copy) %in% samples),]
      breaks <- c(breaks, breaks[length(breaks)] + length(samples))
      cat(breaks, samples,"\n")
    }
    breaks <- c(breaks, nrow(stmatrix))
    ordering <- c(ordering, rownames(stmatrix.copy))
    stmatrix <- stmatrix[ordering,]
  }

  
  mx <- nrow(stmatrix)
  if(is.null(mx)){ 
    stmatrix <- t(t(stmatrix))
    mx <- nrow(stmatrix)
  }
  
  plot(c(1, mx), c(-0.15, 1), t = 'n', xaxt='n', xlab = "", ylab = "Cluster membership (%)", yaxt='n', main=paste0("STRUCTURE, K = ", ncol(stmatrix)))
  dsum <- rep(0, mx)
  mcol <- 2
  apply(stmatrix, 2, function(x){
    for(i in 1:length(x)){
      rect(i-0.5, dsum[i], i+0.5, dsum[i] + x[i], col=brewer.pal(6,"Accent")[mcol], border = "white",lwd=0.1,density=NA)
    }
    dsum <<- dsum + x
    mcol <<- mcol + 1
    cat(dsum, "\n")
  })
  text(x = 1:mx, y = rep(-0.05, mx), popinfo[rownames(stmatrix)], srt = 90, cex=0.8,col=cols[popinfo[rownames(stmatrix)]])
  mids <- diff(breaks) / 2 + 0.5
  for(x in 1:length(mids)) mids[x] <- mids[x] + breaks[x]
  
  last <- mids[length(mids)]
  mids <- mids[-length(mids)]
  
  text(x = mids, y = rep(-0.1, length(mids)-1), paste0("Cluster ", 1:length(mids-1)), cex=1) 
  text(x = last, y = -0.1, paste0(""), cex=1) 
  #axis(1, at=1:nrow(stmatrix), rownames(stmatrix), las = 2, cex.axis = 1)
  axis(2, at=seq(0, 1, 0.1), seq(0, 100, 10), las = 2, cex.axis = 1, lwd=0.2)
  abline(v = breaks + 0.5, lwd=0.2, lty=3)
  box(lwd=0.2)
  return(stmatrix)
  return(breaks)
}

library(StAMPP)
  absnpdata <- read.csv("filtered_snps_AB_NO_DN2.txt", sep="\t", check.names=FALSE)
stammpinput <- t(absnpdata)
stammpinput <- cbind(Sample = rownames(stammpinput), Pop = as.character(samples[rownames(stammpinput),"Breed"]), Ploidy = 2, Format = "BiA", stammpinput)
stammpinput <- as.data.frame(stammpinput)

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values
stammp.D.ind <- stamppNeisD(stammpinput.freq, FALSE) # Population D values

numsnpdata <- read.csv("filtered_snps_numeric_NO_DN2.txt", sep="\t", check.names=FALSE)
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


# Create colors
cols <- c("red", "blue", "orange", "black")
names(cols) <- c("T", "D", "Ni", "Nu")

viewn <- c("T", "D", "Ni", "Nu")
names(viewn) <- c("Tagg", "Dese", "Ni", "Nu")

labelCol <- function(x) {
  if (is.leaf(x)) {
    hclass <- as.character(attr(x, "label"))             # Fetch the class label
    hcol <- cols[hclass]                            # Determine color of the label
    cat(attr(x, "label"), hclass, hcol, "\n")
    attr(x, "nodePar") <- list(lab.col=hcol,labels.cex=2)
  }
  return(x)
}


clustering1$labels <- viewn[clustering1$labels]
clustering2$labels <- viewn[clustering2$labels]

dendrogram1 <- as.dendrogram(clustering1)
dendrogram2 <- as.dendrogram(clustering2)
dendrogram1.col <- dendrapply(dendrogram1, labelCol)
dendrogram2.col <- dendrapply(dendrogram2, labelCol)



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

getMAFs <- function(x){
  return(round(c(length(which(x <= 0.05)) / length(x),
  length(which(x > 0.05 & x <= 0.1)) / length(x),
  length(which(x > 0.1 & x <= 0.3)) / length(x),
  length(which(x > 0.3 & x <= 0.5)) / length(x)),3))
}


png(paste0("Fig1.png"), width = 2200, height = 1200, res=600, pointsize = 3,family = "Arial")
layout(matrix(c(1,1,1,1,2,2,3,4,2,2,5,6), 3, 4, byrow = TRUE), widths=c(1, 0.5, 0.5, 0.5), heights=c(2,1,1))
op <- par(mar=c(5, 4, 4, 2) + 0.1)
op <- par(lwd=0.2)

#op <- par(cex=1.4)
plot(dendrogram2.col, main = "Nei's genetic distance", cex.main=1, cex.axis=1, las=2, ylab="Distance",yaxt='n')
axis(2, seq(0, 0.35, 0.05), seq(0, 0.35, 0.05), las = 2, cex.axis = 1, lwd=0.2)
#op <- par(cex=1)
setwd("D:/Edrive/Goat/DNA/SihamAnalysis/")
## Structure results (run Structure)
results <- paste("D:/Edrive/Goat/DNA/Structure/100k/Goat k",2:5,"/100k/Results/100k_run_1_f",sep="")
op <- par(mar=c(2,4,2,1))
trim <- function (x) gsub("^\\s+|\\s+$", "", x)       # Returns string w/o leading or trailing whitespace
cnt <- 2
for(analysis in results[1]){
  structuredata <- readLines(analysis, warn = FALSE)
  bt <- which(grepl("Estimated Ln Prob of Data", structuredata))                        # Start of info
  st <- which(grepl("Inferred ancestry of individuals:", structuredata))+2              # Start of data table
  et <- which(grepl("Estimated Allele Frequencies in each cluster", structuredata))-3   # End of data table
  datalines <- strsplit(trim(unlist(lapply(strsplit(structuredata[st:et], ":"), "[", 2))), " ")
  stmatrix <- NULL                                              # Structure matrix, holding group % for each group
  for(l in datalines){
    stmatrix <- rbind(stmatrix, as.numeric(strsplit(l, " ")))
  }
  rownames(stmatrix) <- unlist(lapply(strsplit(structuredata[st:et], "\""),"[", 2))
  #png(paste0("Structure100k_K", cnt, ".png"), width = 1400, height = 1000)
    aa <- plotStructure(stmatrix, popinfo, TRUE, FALSE)
  #dev.off()
  cnt <- cnt+1
}

colz <- brewer.pal(4,"Accent")

op <- par(mar=c(3,4,2,1))
h = hist(MAFs[,"Tagg"], breaks=c(0,0.05, 0.1, 0.3, 0.5), plot=FALSE)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=colz, main="Taggar", xlab="", ylab="",yaxt='n', lwd=0.2, cex.axis = 2, cex.main=2)
axis(2, at=seq(0, 50, 10), seq(0, 50, 10), las = 2, cex.axis = 2, lwd=0.2)
legend("topleft", c("Rare", "Intermediate", "Common", "Frequent"), fill =colz, bty='y', bg="white",cex=1.5)

h = hist(MAFs[,"Dese"], breaks=c(0,0.05, 0.1, 0.3, 0.5), plot=FALSE)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=colz, main="Desert", xlab="", ylab="",yaxt='n', lwd=0.2, cex.axis = 2, cex.main=2)
axis(2, at=seq(0, 50, 10), seq(0, 50, 10), las = 2, cex.axis = 2, lwd=0.2)
legend("topleft", c("Rare", "Intermediate", "Common", "Frequent"), fill =colz, bty='y', bg="white",cex=1.5)

h = hist(MAFs[,"Ni"], breaks=c(0,0.05, 0.1, 0.3, 0.5), plot=FALSE)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=colz, main="Nilotic", xlab="", ylab="",yaxt='n', lwd=0.2, cex.axis = 2, cex.main=2)
axis(2, at=seq(0, 50, 10), seq(0, 50, 10), las = 2, cex.axis = 2, lwd=0.2)
legend("topleft", c("Rare", "Intermediate", "Common", "Frequent"), fill = colz, bty='y', bg="white",cex=1.5)

h = hist(MAFs[,"Nu"], breaks=c(0,0.05, 0.1, 0.3, 0.5), plot=FALSE)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, col=colz, main="Nubian", xlab="", ylab="",yaxt='n', lwd=0.2, cex.axis = 2, cex.main=2)
axis(2, at=seq(0, 50, 10), seq(0, 50, 10), las = 2, cex.axis = 2, lwd=0.2)
legend("topleft", c("Rare", "Intermediate", "Common", "Frequent"), fill =colz, bty='y', bg="white",cex=1.5)
dev.off()


op <- par(mfrow=c(2,2))
cnt <- 2
for(analysis in results){
  structuredata <- readLines(analysis, warn = FALSE)
  bt <- which(grepl("Estimated Ln Prob of Data", structuredata))                        # Start of info
  st <- which(grepl("Inferred ancestry of individuals:", structuredata))+2              # Start of data table
  et <- which(grepl("Estimated Allele Frequencies in each cluster", structuredata))-3   # End of data table
  datalines <- strsplit(trim(unlist(lapply(strsplit(structuredata[st:et], ":"), "[", 2))), " ")
  stmatrix <- NULL                                              # Structure matrix, holding group % for each group
  for(l in datalines){
    stmatrix <- rbind(stmatrix, as.numeric(strsplit(l, " ")))
  }
  rownames(stmatrix) <- unlist(lapply(strsplit(structuredata[st:et], "\""),"[", 2))
  #png(paste0("Structure100k_K", cnt, ".png"), width = 1400, height = 1000)
    aa <- plotStructure(stmatrix, popinfo, TRUE, FALSE)
  #dev.off()
  cnt <- cnt+1
}

#hist(MAFs[,"Dese"], breaks=c(0,0.05, 0.1, 0.3, 0.5), col=c("red", "purple", "blue", "green"), main="Desert", xlab="Allele frequency", freq=TRUE)
#legend("topleft", c("Rare", "Intermediate", "Common", "Very common"), fill =c("red", "purple", "blue", "green"), bty='n')
#hist(MAFs[,"Ni"], breaks=c(0,0.05, 0.1, 0.3, 0.5), col=c("red", "purple", "blue", "green"), main="Nilotic", xlab="Allele frequency", freq=TRUE)
#legend("topleft", c("Rare", "Intermediate", "Common", "Very common"), fill =c("red", "purple", "blue", "green"), bty='n')
#hist(MAFs[,"Nu"], breaks=c(0,0.05, 0.1, 0.3, 0.5), col=c("red", "purple", "blue", "green"), main="Nubian", xlab="Allele frequency", freq=TRUE)
#legend("topleft", c("Rare", "Intermediate", "Common", "Very common"), fill =c("red", "purple", "blue", "green"), bty='n')

# Analysis of the PAT and MAT TRD regions in BFMI mice using BxD QTL data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2017
# first written ??

setwd("D:/Edrive/Mouse/DNA/annotation/")
chrinfo <- read.table("chrLengths.txt",sep="\t", header=F, colClasses=c("character", "numeric"))

chrinfo <- chrinfo[-21,] # Ignore the Y chromosome

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
patRegions <- read.table("Analysis/ATB_PAT.txt",sep="\t",header=TRUE, check.names=FALSE, colClasses="character")
matRegions <- read.table("Analysis/ATB_MAT.txt",sep="\t",header=TRUE, check.names=FALSE, colClasses="character")

allRegions <- rbind(patRegions, matRegions)

setwd("D:/Edrive/Mouse/DNA/MegaMuga/GeneNetwork")
searchres <- read.csv("BXDpheno.txt", sep = "\t", header=FALSE, colClasses="character")

chr <- unlist(lapply(strsplit(searchres[,7], ": "),"[", 1))
loc <- as.numeric(as.character(unlist(lapply(strsplit(searchres[,7], ": "),"[", 2)))) * 1000000

chr <- gsub("Chr", "", chr)

searchres <- cbind(searchres, chr=chr, loc=loc)
searchres <- searchres[-which(searchres[,"chr"] == "N/A"),]

searchres[1:5,]

significant.lrs <- which(as.numeric(searchres[, 6]) > 16) # LRS scores not LOD :(
searchres <- searchres[significant.lrs,]

testQTLs <- function(allRegions, keywords = c("Obesity", "Weight")){
  nQTL <- 0
  nKeyw <- 0
  for(x in 1:nrow(allRegions)){
    idx <- which(searchres[, "chr"] == allRegions[x,"Chr"] & as.numeric(searchres[, "loc"]) > (as.numeric(allRegions[x,"Start"])- 500000) & searchres[, "loc"] < (as.numeric(allRegions[x,"Stop"])+ 500000))
    hasKW <- FALSE
    hasQTL <- FALSE
    if(length(idx) > 0){
      hasQTL <- TRUE
      for(word in keywords){
        grepidx <- which(grepl(word, searchres[idx,2], ignore.case=TRUE))
        if(length(grepidx) > 0) {
          cat(x, " - ", word, " - ", idx[grepidx],"\n")
          for(ii in idx[grepidx]){
            cat(unlist(allRegions[x,]), unlist(searchres[ii,]),"\n",file="regionsQTLMat.txt",append=TRUE,sep="\t")
          }
          hasKW <- TRUE
        }
      }
    }
    if(hasQTL) nQTL <- nQTL + 1
    if(hasKW) nKeyw <- nKeyw + 1
  }
  cat("Regions with a QTL:", nQTL / nrow(allRegions) * 100, "\n")
  cat("Keywords per QTL:", nKeyw / nQTL * 100, "\n")
  cat("Keywords per region:", nKeyw / nrow(allRegions) * 100, "\n")
  return(list(nQTL, nKeyw))
}

randomRegions <- function(allRegions, chrinfo){
  newRegions <- allRegions
  for(x in 1:nrow(allRegions)){
    newChr <- sample(chrinfo[,1])[1]
    chrLength <- as.numeric(chrinfo[chrinfo[,1] == newChr,2])
    rL <- as.numeric(as.character(allRegions[x,"Stop"])) - as.numeric(as.character(allRegions[x,"Start"]))
    newStart <-  round(runif(1, 0, chrLength - rL))
    newStop <- newStart + rL
    newRegions[x, "Chr"] <- newChr
    newRegions[x, "Start"] <- newStart
    newRegions[x, "Stop"] <- newStop
  }
  return(newRegions)
}

realpat <- testQTLs(patRegions)
realmat <- testQTLs(matRegions)
realCombi <- testQTLs(allRegions)

nperms <- 10000

pPat <- vector("list", nperms)
pMat <- vector("list", nperms)
pCombi <- vector("list", nperms)
x <- 1
for(x in 1:nperms){
  newRegions <- randomRegions(patRegions,chrinfo)
  pPat[[x]] <- testQTLs(newRegions)

  newRegions <- randomRegions(matRegions,chrinfo)
  pMat[[x]] <- testQTLs(newRegions)

  newRegions <- randomRegions(allRegions,chrinfo)
  pCombi[[x]] <- testQTLs(newRegions)
}

# Mean and sd for n QTLs in 85 regions
cat(mean(unlist(lapply(pPat,"[",1))/nrow(patRegions)), sd(unlist(lapply(pPat,"[",1))/nrow(patRegions)), "\n")
cat(mean(unlist(lapply(pMat,"[",1))/nrow(matRegions)), sd(unlist(lapply(pMat,"[",1))/nrow(matRegions)), "\n")
cat(mean(unlist(lapply(pCombi,"[",1))/nrow(allRegions)), sd(unlist(lapply(pCombi,"[",1))/nrow(allRegions)), "\n")

pvalQTLsPAT <- (length(which(realpat[[1]][1] < sort(unlist(lapply(pPat,"[",1))))) + 1) / nperms
pvalQTLsMAT <- (length(which(realmat[[1]][1] < sort(unlist(lapply(pMat,"[",1))))) + 1) / nperms
pvalQTLsCombi <- (length(which(realCombi[[1]][1] < sort(unlist(lapply(pCombi,"[",1))))) + 1) / nperms

cat(pvalQTLsPAT, pvalQTLsMAT, pvalQTLsCombi, "\n")

# Mean and sd for n Keywords in 85 regions
cat(mean(unlist(lapply(pPat,"[",2))/nrow(patRegions)), sd(unlist(lapply(pPat,"[",2))/nrow(patRegions)), "\n")
cat(mean(unlist(lapply(pMat,"[",2))/nrow(matRegions)), sd(unlist(lapply(pMat,"[",2))/nrow(matRegions)), "\n")
cat(mean(unlist(lapply(pCombi,"[",2))/nrow(allRegions)), sd(unlist(lapply(pCombi,"[",2))/nrow(allRegions)), "\n")

pvalKWsPAT <- (length(which(realpat[[2]][1] < sort(unlist(lapply(pPat,"[",2))))) + 1) / nperms
pvalKWsMAT <- (length(which(realmat[[2]][1] < sort(unlist(lapply(pMat,"[",2))))) + 1) / nperms
pvalKWsCombi <- (length(which(realCombi[[2]][1] < sort(unlist(lapply(pCombi,"[",2))))) + 1) / nperms

cat(pvalKWsPAT, pvalKWsMAT, pvalKWsCombi, "\n")

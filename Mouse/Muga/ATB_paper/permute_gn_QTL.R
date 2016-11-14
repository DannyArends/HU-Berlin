setwd("E:/Mouse/DNA/MegaMuga/")

patRegions <- read.table("Analysis/ATB_PAT.txt",sep="\t",header=TRUE, check.names=FALSE, colClasses="character")
matRegions <- read.table("Analysis/ATB_MAT.txt",sep="\t",header=TRUE, check.names=FALSE, colClasses="character")

#allRegions <- matRegions
allRegions <- rbind(matRegions, patRegions)

setwd("E:/Mouse/DNA/MegaMuga/GeneNetwork")
searchres <- read.csv("BXDpheno.txt", sep = "\t", header=FALSE, colClasses="character")

chr <- unlist(lapply(strsplit(searchres[,7], ": "),"[", 1))
loc <- as.numeric(as.character(unlist(lapply(strsplit(searchres[,7], ": "),"[", 2)))) * 1000000

chr <- gsub("Chr", "", chr)

searchres <- cbind(searchres, chr=chr, loc=loc)
searchres <- searchres[-which(searchres[,"chr"] == "N/A"),]

searchres[1:5,]

keywords <- c("Obesity", "Metabolism", "Body Weight", "Weight", "Fat", "Growth")

chrs <- as.character(unique(searchres[,"chr"]))
res <- NULL
for(x in chrs){
  res <- rbind(res, c(x, max(as.numeric(as.character(searchres[searchres[,"chr"] == x, "loc"])))))
}

testQTLs <- function(allRegions){
  nQTL <- 0
  nKeyw <- 0
  for(x in 1:nrow(allRegions)){
    idx <- which(searchres[, "chr"] == allRegions[x,"Chr"] & as.numeric(searchres[, "loc"]) > (as.numeric(allRegions[x,"Start"])- 500000) & searchres[, "loc"] < (as.numeric(allRegions[x,"Stop"])+ 500000))
    hasKW <- FALSE
    hasQTL <- FALSE
    if(length(idx) > 0){
      hasQTL <- TRUE
      for(word in keywords){
        if(any(grepl(word, searchres[idx,2], ignore.case=TRUE))){
          #cat(" - ", word,"\n")
          hasKW <- TRUE
        }
      }
    }
    if(hasQTL) nQTL <- nQTL + 1
    if(hasKW) nKeyw <- nKeyw + 1
  }
  cat("Regions with QTL:", nQTL / nrow(allRegions) * 100, "\n")
  cat("Keywords per QTL:", nKeyw / nQTL * 100, "\n")
  return(list(nQTL, nKeyw))
}

randomRegions <- function(allRegions){
  newRegions <- allRegions
  for(x in 1:nrow(allRegions)){
    newChr <- sample(chrs)[1]
    rL <- as.numeric(as.character(allRegions[x,"Stop"])) - as.numeric(as.character(allRegions[x,"Start"]))
    newStart <-  round(runif(1, 0, rL))
    newStop <- newStart + rL
    newRegions[x, "Chr"] <- newChr
    newRegions[x, "Start"] <- newStart
    newRegions[x, "Stop"] <- newStop
  }
  return(newRegions)
}

testQTLs(allRegions)
res <- vector("list", 1000)
x <- 1
for(x in 1:1000){
  newRegions <- randomRegions(allRegions)
  res[[x]] <- testQTLs(newRegions)
}

mean(unlist(lapply(res,"[",2)) / unlist(lapply(res,"[",1)))
sd(unlist(lapply(res,"[",2)) / unlist(lapply(res,"[",1)))
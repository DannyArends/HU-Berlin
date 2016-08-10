# Analysis of STRUCTURE results
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jul, 2016
# first written Jul, 2016

setwd("E:/Goat/DNA/SihamAnalysis")
samples <- read.table("sampleinfo.txt", sep="\t")
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`
samples    <- cbind(samples, locations[rownames(samples),])
samples    <- cbind(samples, locationShort = as.character(unlist(lapply(strsplit(as.character(samples[,"Location"]), "_"), "[",1))))
samples    <- samples[-which(rownames(samples) == "DN 2"),] # Throw away the duplicate individual because it confuses STRUCTURE

popinfo <- as.character(samples[,"Breed"])
names(popinfo) <- gsub(" ", "", rownames(samples))


##E:\Goat\DNA\Structure\test\10k2.5k\Results
## Structure results (run Structure)
structuredir <- "E:/Goat/DNA/Structure/"
projectname <- "test"
paramsetname <- "10k2.5k"

fullpath <- paste0(structuredir, "/", projectname, "/", paramsetname, "/Results")
structureruns <- dir(fullpath)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)       # Returns string w/o leading or trailing whitespace

analyzeStructure <- function(stmatrix, popinfo, confidence = 0.0){
  stmatrix <- stmatrix[which(apply(stmatrix,1,max) > confidence), ]
  mp <- apply(stmatrix, 1, which.max)
  mp <- cbind(mp, Strain = NA)

  ii <- rownames(mp)[which(rownames(mp) %in% names(popinfo))]
  mp[ii, "Strain"] <- as.character(popinfo[ii])
  mp[is.na(mp)] <- "X"

  counts <- matrix(0, ncol(stmatrix), 4, dimnames=list(1:ncol(stmatrix), c("H","S","K", "P")))
  for(x in 1:nrow(mp)){
    counts[as.numeric(mp[x,1]), mp[x,2]] <- counts[mp[x,1], mp[x,2]] + 1
  }
  print(counts)
  return(counts)
}

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
  
  plot(c(1, mx), c(-0.15, 1), t = 'n', xaxt='n', xlab = "Individual", ylab = "Cluster membership (%)", yaxt='n', main=paste0("STRUCTURE, clusters = ", ncol(stmatrix)))
  dsum <- rep(0, mx)
  mcol <- 2
  apply(stmatrix, 2, function(x){
    for(i in 1:length(x)){
      rect(i-0.4, dsum[i], i+0.4, dsum[i] + x[i], col=mcol)
    }
    dsum <<- dsum + x
    mcol <<- mcol + 1
    cat(dsum, "\n")
  })
  text(x = 1:mx, y = rep(-0.05, mx), popinfo[rownames(stmatrix)], srt = 90)
  mids <- diff(breaks) / 2 + 0.5
  for(x in 1:length(mids)) mids[x] <- mids[x] + breaks[x]
  
  last <- mids[length(mids)]
  mids <- mids[-length(mids)]
  
  text(x = mids, y = rep(-0.1, length(mids)-1), paste0("Cluster ", 1:length(mids-1)), cex=0.8) 
  text(x = last, y = -0.1, paste0("Unassigned"), cex=0.8) 
  axis(1, at=1:nrow(stmatrix), rownames(stmatrix), las = 2, cex.axis = 0.8)
  axis(2, at=seq(0, 1, 0.1), seq(0, 100, 10), las = 2, cex.axis = 0.8)
  abline(v = breaks + 0.5,lwd=0.5, lty=2)
  return(stmatrix)
  return(breaks)
}

cnt <- 1
for(analysis in paste0(fullpath, "/", structureruns)){
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
  cat("--", analysis, "--\n")
  cat(structuredata[bt:(st-5)],sep="\n")        # Print some structure information
  png(paste0(analysis, ".st.png"), width = 1024, height = 800)
  aa <- plotStructure(stmatrix, popinfo, TRUE, FALSE)
  dev.off()
  #cat("--All individuals--\n")
  #counts <- analyzeStructure(stmatrix)          # Compare how good the structure model fits with the breeders perspective
  #cat("--High purity (> 0.8)--\n")
  #counts <- analyzeStructure(stmatrix, 0.8)     # Compare how good the pure structure model fits with the breeders perspective
  #cat("----\n")
  #if(ncol(stmatrix) == 4) break                 # We did more but I want too see if Saria's hypothesis is correct
  cnt <- cnt +1
}

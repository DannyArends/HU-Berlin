#
# Server script to extract SNPs for a region
#

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); } #q("no") }
}

samples <- c("ext_L7254", "ext_L7256", "ext_L7257", "ext_L7255", "ext_L7258")
names(samples) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")

## Start of analysis

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)

#setwd("~/DIFE/")
setwd("E:/Mouse/DNA/Sequencing/DifeMouse/")
regions <- read.table("regions_bfmi.txt", sep = "\t",header=TRUE, colClasses="character", row.names=1)   # Start loading the BFMi regions

QTLregions <- regions[which(regions[,"type"] == "QTL"),]
DivAregions <- regions[which(regions[,"type"] == "DivA"),]
TRDregions <- regions[which(regions[,"type"] == "TRD"),]

# test if 2 is within 1
within <- function(s1, s2, e1, e2, wiggle = 500000) {
  if(s2 >= s1) {
    if(e2 <= (e1 + wiggle)) return(TRUE);
  }
  if(e2 <= e1) {
    if(s2 >= (s1 - wiggle)) return(TRUE);
  }
  return(FALSE)
}

# Analysis of QTL regions

getInsideMatrix <- function(chrData) {
  if(nrow(chrData) > 1){
    insides <- matrix(0, nrow(chrData),nrow(chrData))
    rownames(insides) <- rownames(chrData)
    colnames(insides) <- rownames(chrData)
    for(x in 1:nrow(chrData)){
      for(y in 1:nrow(chrData)) {
        inside <- within(as.numeric(chrData[x,"start"]), 
                   as.numeric(chrData[y,"start"]), 
                   as.numeric(chrData[x,"end"]), 
                   as.numeric(chrData[y,"end"]))
        insides[x,y] <- as.numeric(inside)
      }
    }
    return(insides)
  }
  return(matrix(0, 1,1))
}

minimumregions <- NULL
for(chr in unique(QTLregions[,"chr"])){
  chrData <- QTLregions[which(QTLregions[,"chr"] == chr),]
  insideMatrix <- getInsideMatrix(chrData)
  print(insideMatrix)
  if(nrow(insideMatrix) > 0){
    cat("CHR", chr, "\n")
    for(x in 1:nrow(insideMatrix)){
      idx <- which(insideMatrix[x,] == 1)
      cat(rownames(insideMatrix)[x], " ")
      if(length(idx) > 0){
        insideX <- rownames(insideMatrix)[idx]
        lengths <- as.numeric(chrData[insideX, "end"]) - as.numeric(chrData[insideX, "start"])
        cat("min:", insideX[which.min(lengths)], "\n")
        if(insideX[which.min(lengths)] != rownames(insideMatrix)[x]){         # if the miniumum region is not the current region
          minimumregions <- c(minimumregions, insideX[which.min(lengths)])
        }
      }
      cat(x, idx,"\n")
    }
  }
}

mcolors <- c("142", "144", "purple", "146", "148", "150", "502", "102", "258", "635", "637", "420", "27", "45", "53")
names(mcolors) <- c("FAT", "FAT%", "LEAN", "Gonadal AT", "Renal AT", "Inguinal AT", "BAT", "Liver", "BW", "RepF", "SubF", "BG", "Insulin", "TG", "Chol")

smallestSubSet <- unique(minimumregions)
smallesRegionSet <- QTLregions[smallestSubSet,]
smallesRegionSet <- smallesRegionSet[sort(as.numeric(smallesRegionSet[,"chr"]), index.return=T)$ix,]

smallesRegionSet <- cbind(smallesRegionSet, RegionsOverlap = NA)

op <- par(mfrow = c(2, (nrow(smallesRegionSet)+1)/2))
op <- par(mar = c(5,3,2,2))
for(x in 1:nrow(smallesRegionSet)){
  chr <- smallesRegionSet[x, "chr"]
  chrData <- QTLregions[which(QTLregions[,"chr"] == chr),]
  insideMatrix <- getInsideMatrix(chrData)
  sName <- rownames(smallesRegionSet)[x]
  biggerThenS <- rownames(insideMatrix)[which(insideMatrix[,sName] == 1)]
  lengthBig <- round((as.numeric(chrData[biggerThenS, "end"]) - as.numeric(chrData[biggerThenS, "start"])) / 1000000,1)
  lengthS <- round((as.numeric(chrData[sName, "end"]) - as.numeric(chrData[sName, "start"])) / 1000000,1)
  smallesRegionSet[x,"RegionsOverlap"] <- paste0(biggerThenS, collapse = ", ")
  ordering <- sort(lengthBig, index.return = TRUE, decreasing = TRUE)$ix
  cat(paste0(biggerThenS,"-",lengthBig)[ordering], "->>", sName,"-",lengthS, "\n")
  chrL <- chrInfo[chrInfo[,"Chr"] == chr,"Length"]
  
  
  plot(c(0.75,length(ordering)+1.25), c(0, chrL), t = 'n', xlab="", ylab="position", yaxt='n', xaxt='n', main=paste0("Chr",chr))
  cnt <- 1
  for(rName in biggerThenS[ordering]){
    ys <- as.numeric(chrData[rName, "start"])
    ye <- as.numeric(chrData[rName, "end"])
    points(c(cnt,cnt), c(ys, ye), t = 'l', lwd = 2, col=mcolors[chrData[rName, "Phenotype"]])
    cnt <- cnt + 1
  }
  ys <- as.numeric(chrData[sName, "start"])
  ye <- as.numeric(chrData[sName, "end"])
  points(c(cnt,cnt), c(ys, ye), t = 'l', lwd = 2, col=mcolors[chrData[rName, "Phenotype"]])
  
  axis(2, at = seq(0, chrL, 10000000), seq(0, chrL, 10000000) / 1000000, las=2)
  axis(1, at = 1:(length(ordering)+1), c(biggerThenS,sName), las=2, cex.axis=0.8)
}

plot(c(0,2), c(1,16), t = 'n', yaxt='n', xaxt='n', xlab="", ylab="")
for(y in 1:16){
  points(c(0,0.5), c(y,y), col=mcolors[y], lwd=2, t='l')
  text(1.5, y, names(mcolors)[y])
}

write.table(smallesRegionSet, "smallestregionsubset.txt", sep="\t", quote = FALSE)
 
# SNP analysis
 
# Use tabix to subset the big vcf.gz file
for(i in 1:nrow(regions)){
  outfile <- paste0("~/DIFE/analysis/", regions[i,"id"], ".txt")
  cmd <- paste0("tabix -h /home/arends/NAS/Mouse/DNA/Sequencing/DifeMouse/RAW/ALL_variants.vcf.gz chr", regions[i,"chr"], ":", regions[i,"start"], "-", regions[i,"end"], " > ", outfile)
  execute(cmd, outfile)
}


### OLD code from here ###

# Plot all regions for a certain chromosome
plotRegions <- function(regions, chr = 1){
  chrsubset <- regions[which(regions[,"chr"] == chr),]
  mY <- max(chrsubset[, "end"], na.rm=TRUE)
  mX <- nrow(chrsubset)
  plot(c(0,mX), c(0,mY), t = 'n', yaxt='n')
  for(i in 1:nrow(chrsubset)){
    points(c(i,i), c(chrsubset[i,"start"], chrsubset[i, "end"]), t = 'l', col=i,lwd=2)
  }
  axis(2, seq(0, mY, 10000000), seq(0, mY, 10000000)/ 1000000, las=2)
}

plotRegions(regions, 3)

loadRegion <- function(name = "HUregion1", loc = "~/DIFE/analysis/") {
  header <- strsplit(readLines(paste0(loc, name, ".txt"), n = 133)[133], "\t")[[1]]
  header[1] <- "CHROM"
  regiondata <- read.table(paste0(loc, name, ".txt"),sep="\t")
  colnames(regiondata) <- header
  idx <- which(colnames(regiondata) %in% samples)
  colnames(regiondata)[idx] <- names(samples)[match(colnames(regiondata)[idx], samples)]

  snpnames <- as.character(regiondata[,"ID"])
  noID <- which(snpnames == ".")
  snpnames[noID] <- paste0("N-", gsub("chr", "", regiondata[noID, "CHROM"]), ":", regiondata[noID, "POS"])
  rownames(regiondata) <- snpnames
  return(regiondata)
}

onlyGeno <- function(regiondata){
  geno <- t(apply(regiondata[,names(samples)], 1, function(x){
    unlist(lapply(strsplit(x, ":"), "[", 1))
  }))
  rownames(geno) <- rownames(regiondata)
  return(geno)
}

filterGeno <- function(genotypes, mfilter, naIsOK = FALSE){
  dontcare <- names(which(is.na(mfilter)))
  geno <- genotypes[, -which(colnames(genotypes) %in% dontcare)]
  mfilter <- mfilter[-which(is.na(mfilter))]
  
  matches <- NULL
  for(x in 1:nrow(geno)){
    if(all(geno[x, names(mfilter)] == mfilter)){
      matches <- cbind(matches, x)
    }
  }
  return(genotypes[matches,])
}

regiondata <- loadRegion()
genotypes <- onlyGeno(regiondata)

mfilterS1 <- c(NA, "1/1", "0/0", NA, NA); names(mfilterS1) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")
mfilterS2 <- c(NA, "0/0", "1/1", NA, NA); names(mfilterS2) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")

S1S2diff <- rbind(filterGeno(genotypes, mfilterS1), filterGeno(genotypes, mfilterS2))


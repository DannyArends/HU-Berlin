#
# Analysis of BxD genotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

setwd("E:/Mouse/BxD")
alldata <- read.csv("genotypes_26-05.txt", sep = "\t", skip = 1, row.names=2, header=TRUE, na.strings=c("NA", "-3", ""), colClasses = "character")

map <- alldata[,c(3, 5, 7)]                                                     # extract the physical map
map[,2] <- as.numeric(map[,2])
genotypes <- alldata[,-c(1:11)][,1:203]                                         # extract the genotypes

wrongEncoding <- names(which(lapply(apply(genotypes, 2, table), length) > 3))
wE <- apply(genotypes[,wrongEncoding],2, function(x){ rownames(genotypes)[which(!(x %in% c("-1","0","1")))] })

for(x in 1:length(wE)){
  genotypes[wE[[x]], names(wE)[x]] <- NA
}

toNumGeno <- function(genotypes){
  numgeno <- apply(genotypes, 2, function(x){ return(as.numeric(as.character(x))) })      # transform genotypes to numeric values
  rownames(numgeno) <- rownames(genotypes)                                                # Set the markernames
  return(numgeno)
}

geno <- toNumGeno(genotypes)

chr <- 6
smap <- map[which(map[,"Chr"] == chr), ]
sgeno <- geno[rownames(smap),]
recM <- NULL
for(x in rownames(smap)){
  sLoc <- smap[x, 2] - 0.5 
  eLoc <- smap[x, 2] + 0.5
  ii <- which(smap[,2] > sLoc & smap[,2] < eLoc)
  if(length(ii) > 1){
    recs <- NULL
    for(ind in colnames(sgeno)){
      rec <- 0
      for(index in 1:(length(ii)-1)){
        if(!is.na(sgeno[ii[index], ind]) && !is.na(sgeno[ii[index+1], ind]) && sgeno[ii[index], ind] != sgeno[ii[index+1], ind]){  # Recombination
          rec <- rec + 1
        }
      }
      #if(rec > 1) cat(x, ind, rec, unlist(sgeno[ii, ind]),"\n")
      recs <- c(recs, rec)
    }
    recM <- rbind(recM, recs)
  }
  cat(x, "\n")
}

image(1:nrow(recM), 1:ncol(recM), recM)
axis(2, at=1:ncol(recM), colnames(geno), las=2, cex.axis=0.3)

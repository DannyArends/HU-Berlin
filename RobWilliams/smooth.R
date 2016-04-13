#
# Smoothing genotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

library("openxlsx")
source("~/BxD/functions.R")

setwd("~/BxD/")
alldata <- read.csv("genotypes.txt", sep = "\t", colClasses="character", skip = 1, row.names=2, header=TRUE, na.strings=c("NA", ""))

map <- alldata[,c(3, 5, 7)]                                                       # Extract the physical map
genotypes <- alldata[,-c(1:11)][,1:203]                                           # Extract the genotypes

genotypes[3238, "BXD53"]  <- NA
genotypes[3239, "BXD53"]  <- NA
genotypes[3240, "BXD53"]  <- NA
genotypes[3241, "BXD53"]  <- NA
genotypes[20364,"BXD146"] <- NA
genotypes[16041,"BXD93"]  <- NA
genotypes[16042,"BXD93"]  <- NA

#table(unlist(genotypes))

toNumGeno <- function(genotypes){
  numgeno <- apply(genotypes, 2, function(x){ return(as.numeric(as.character(x))) })      # transform genotypes to numeric values
  rownames(numgeno) <- rownames(genotypes)                                                # Set the markernames
  return(numgeno)
}

geno <- toNumGeno(genotypes)

geno <- geno[-which(rownames(geno) == "Affy_10542764_Arntl2"),]             # Marker without a locations
map <- map[-which(rownames(map) == "Affy_10542764_Arntl2"),]                # Marker without a locations

#table(unlist(geno))

recomb.index <- apply(geno, 2, function(x){recombinations(x, chr = map[,1]) })
recomb.loc <- apply(geno, 2, function(x){recombinations(x, loc = as.numeric(map[,2]), chr = map[,1], TRUE) })

recomb.changed <- NULL
ngeno <- geno
for(i in 1:ncol(geno)){
  diffs <- diff(as.numeric(get.locs(recomb.loc[[i]])))
  chrs  <- diff(as.numeric(as.factor(get.chr(recomb.loc[[i]]))))
  doubleR <- which(diffs > 0 & diffs < 1 & chrs == 0)
  for(x in doubleR){
    markerD <- recomb.index[[i]][x + 1] - recomb.index[[i]][x]  
    s <- floor(recomb.index[[i]][x])
    while(is.na(geno[s,i])){ s <- s - 1; }
    sb <- s
    e <- ceiling(recomb.index[[i]][x+1])
    while(is.na(geno[e,i])){ e <- e + 1; }
    eb <- e
    while(!is.na(geno[sb-1,i]) && sb > 1 && map[s, "Chr"] == map[sb-1, "Chr"] && geno[s,i] == geno[sb-1,i]){
      sb <- sb-1
    }
    while(!is.na(geno[eb+1,i]) && eb < nrow(geno) && map[e, "Chr"] == map[eb+1, "Chr"] && geno[e,i] == geno[eb+1,i]){
      eb <- eb+1
    }
    if(geno[s,i] == geno[e,i]){  # Double recombination between similar XXXXX YYY XXXXX
      if((s-sb) >= 10 && (eb-e) >= 10){
        ngeno[(s+1):(e-1), i] <- rep(geno[s,i], length((s+1):(e-1)))
        recomb.changed <- rbind(recomb.changed, c(i, (s+1), (e-1)))
      }
    }
    if(geno[s,i] != geno[e,i]){  # Double recombination between dissimilar XXXXX YYY ZZZZZ
      if((s-sb) >= 10 && (eb-e) >= 10 && e-s < 10){
        ngeno[(s+1):(e-1), i] <- NA
        recomb.changed <- rbind(recomb.changed, c(i, (s+1), (e-1)))
      }
    }
    cat("Double recombination of length", markerD, s, s-sb, e, eb-e,"|", unlist(geno[s:e,i]),"->",unlist(ngeno[s:e,i]), "\n")
  }
}

geno <- ngeno

recomb.locAfter <- apply(geno, 2, function(x){recombinations(x, loc = as.numeric(map[,2]), chr = map[,1], TRUE) })

nrecombs.before <- lapply(recomb.loc, length)
nrecombs.after <- lapply(recomb.locAfter, length)

write.table(cbind(nrecombs.before, nrecombs.after), "recoms.txt",sep="\t")

wb <- createWorkbook()
sheet <- addWorksheet(wb, sheetName="Genotypes BXD")
writeData(wb, sheet = 1, rownames(geno), startCol = 1, startRow = 2)
writeData(wb, sheet = 1, map[rownames(geno),1], startCol = 2, startRow = 2)
writeData(wb, sheet = 1, map[rownames(geno),2], startCol = 3, startRow = 2)
writeData(wb, sheet = 1, map[rownames(geno),3], startCol = 4, startRow = 2)
for(x in 1:ncol(geno)){ # ncol(geno)
  writeData(wb, sheet = 1, c(colnames(geno)[x], geno[,x]), startCol = (x+4))
}

recombOK <- createStyle(fontColour = rgb(0,0,0), fgFill = rgb(0,0,1))
recombChanged <- createStyle(fontColour = rgb(0,0,0), fgFill = rgb(0.3671096,0.2657807,0.3671096))
recombBAD <- createStyle(fontColour = rgb(0,0,0), fgFill = rgb(1,0,0))

for(i in 1:ncol(geno)){ # ncol(geno)
  rows <- as.numeric(unlist(apply(cbind(floor(as.numeric(recomb.index[[i]])), ceiling(as.numeric(recomb.index[[i]]))), 1, function(x){ return(x[1]:x[2])})))
  addStyle(wb, sheet = 1, recombOK, cols = rep((i+4), length(rows)), rows = (rows+1))
  diffs <- diff(as.numeric(get.locs(recomb.loc[[i]])))
  chrs  <- diff(as.numeric(as.factor(get.chr(recomb.loc[[i]]))))
  BADlocsFrom <- recomb.index[[i]][which(diffs > 0 & diffs < 1 & chrs == 0)]
  BADlocsTo <- recomb.index[[i]][which(diffs > 0 & diffs < 1 & chrs == 0) + 1]
  rows <- as.numeric(unlist(apply(cbind(floor(as.numeric(BADlocsFrom)), ceiling(as.numeric(BADlocsTo))), 1, function(x){ return(x[1]:x[2])})))
  addStyle(wb, sheet = 1, recombBAD, cols = rep((i+4), length(rows)), rows = (rows+1))

  for(x in which(recomb.changed[,1] == i)){
      addStyle(wb, sheet = 1, recombChanged, cols = rep((i+4), length(recomb.changed[x,2]:recomb.changed[x,3])), rows = ((recomb.changed[x,2]:recomb.changed[x,3])+1))
  } 
}

saveWorkbook(wb, file="genotypes_recombinations_2.xlsx", overwrite = TRUE)




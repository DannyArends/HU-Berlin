#
# Smoothing genotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

library("openxlsx")
source("D:/Github/HU-Berlin/RobWilliams/functions.R")

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

table(unlist(geno))
geno[geno == 0] <- NA          # Replace the 0s by missing values (since we don't know the genotype there)
noLoc <- rownames(map)[which(is.na(map[,2]))]

geno <- geno[-which(rownames(geno) %in% noLoc), -c(1:4)]  # remove the no location marker, and the first 4 individuals (founders / F1 individuals)
map <- map[rownames(geno), ]

recomb.index <- apply(geno, 2, function(x){recombinations(x, chr = map[,1]) })
recomb.loc <- apply(geno, 2, function(x){recombinations(x, loc = as.numeric(map[,2]), chr = map[,1], TRUE) })

recomb.changed <- NULL
ngeno <- geno
for(i in 1:ncol(geno)){
  diffs <- diff(as.numeric(get.locs(recomb.loc[[i]])))
  chrs  <- diff(as.numeric(as.factor(get.chr(recomb.loc[[i]]))))
  doubleR <- which(diffs > 0 & diffs < 1.5 & chrs == 0)
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
    naM <- length(which(is.na(ngeno[(s+1):(e-1), i])))  # Number of NA markers in the interval
    if(geno[s,i] == geno[e,i]){  # Double recombination between similar XXXXX YYY XXXXX
      #if((s-sb) >= 5 && (eb-e) >= 5){ # Demand at least 10 before and 10 after that are the same
        
        ngeno[(s+1):(e-1), i] <- rep(geno[s,i], length((s+1):(e-1)))
        recomb.changed <- rbind(recomb.changed, c(i, (s+1), (e-1)))
      #}
    }
    if(geno[s,i] != geno[e,i]){  # Double recombination between dissimilar XXXXX YYY ZZZZZ
      cat("XXXXX YYY ZZZZZ\n")
      if((s-sb) >= 5 && (eb-e) >= 5 && e-s < 10){
        ngeno[(s+1):(e-1), i] <- NA
        recomb.changed <- rbind(recomb.changed, c(i, (s+1), (e-1)))
      }
    }
    cat("[",i,"]Double recombination of length", markerD, naM, s, s-sb, e, eb-e,"|", unlist(geno[s:e,i]),"->",unlist(ngeno[s:e,i]), "\n")
  }
}

recomb.index.n <- apply(ngeno, 2, function(x){recombinations(x, chr = map[,1]) })

plot(unlist(lapply(recomb.loc, length)), pch="+", xaxt = 'n')
points(unlist(lapply(recomb.loc.n, length)), col='green',pch="-",cex=2)
axis(1, 1:ncol(geno), colnames(geno),las=2, cex=0.3)


hist(unlist(lapply(recomb.loc.n, length)), col=rgb(0,1,0,0.6), breaks=seq(0, 200, 10))
hist(unlist(lapply(recomb.loc, length)), col=rgb(0,0,0,0.6), breaks=seq(0, 200, 10), add=TRUE)
hist(unlist(lapply(recomb.loc.n, length)), col=rgb(0,1,0,0.6), breaks=seq(0, 200, 10), add=TRUE)

missingdata <- which(apply(ngeno,1,function(x){length(which(is.na(x)))})/ncol(geno) > 0.1)

ngeno <- ngeno[-missingdata,]
map <- map[-missingdata,]
recomb.index.nf <- apply(ngeno, 2, function(x){recombinations(x, chr = map[,1]) })
recomb.loc.nf <- apply(ngeno, 2, function(x){recombinations(x, loc = as.numeric(map[,2]), chr = map[,1], TRUE) })

plot(c(0,200), c(0,200), t = 'n', xaxt = 'n', xlim=c(0,200),xaxs="i", main="BxD recombinations (0 = Missing)", ylab="nRecombinations", xlab="")
points(unlist(lapply(recomb.loc, length)), pch="+")
points(unlist(lapply(recomb.loc.n, length)), col='red',pch="-",cex=2)
points(unlist(lapply(recomb.loc.nf, length)), col='blue',pch="-",cex=2)

axis(1, 1:ncol(geno), colnames(geno),las=2, cex.axis=0.7)
legend("topright", c("Before ", "After smoothing (1.5 mb window)", "After removal of 'poor quality' markers"), pch=c("+", "-", "-"), lwd=c(1,2,2), col=c("black", "red", "blue"))

mean(unlist(lapply(recomb.loc, length)))      # 64.77387
sd(unlist(lapply(recomb.loc, length)))        # 33.41058

sd(unlist(lapply(recomb.loc.n, length)))      # 27.76904
mean(unlist(lapply(recomb.loc.n, length)))    # 59.35678

mean(unlist(lapply(recomb.loc.nf, length)))   # 54.01508
sd(unlist(lapply(recomb.loc.nf, length)))     # 23.13388

write.table(cbind(alldata[rownames(ngeno), 1:15], ngeno), file="genotypes_27-05_NAs.txt", sep="\t", quote=FALSE, na ="")

ngenoO <- ngeno
ngenoO[is.na(ngenoO)] <- 0

write.table(cbind(alldata[rownames(ngeno), 1:15], ngenoO), file="genotypes_27-05_0s.txt", sep="\t", quote=FALSE, na ="")



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




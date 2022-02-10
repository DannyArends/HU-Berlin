library(strawr)

setwd("C:/Users/Arends/Downloads")
bed <- as.matrix(straw("NONE", "FKR02_1bfmi_L001.allValidPairs.hic", "3", "3", unit = "BP", 100000))
mm <- matrix(0, length(unique(bed[,1])) , length(unique(bed[,2])), dimnames = list(unique(bed[,1]), unique(bed[,2])))
for(x in 1:nrow(bed)){
  mm[as.character(bed[x,1]), as.character(bed[x,2])] <- bed[x,3]
  mm[as.character(bed[x,2]), as.character(bed[x,1])] <- bed[x,3]
}

iix <- as.character(seq(30000000, 40000000, 100000))

image(mm[iix,iix], breaks = c(0, 100, 1000, 1000000), col = c("white", "gray", "black"))

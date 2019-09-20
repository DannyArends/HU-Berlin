#
# Plot results from a Hi-C experiment
#

setwd("D:/")

mold <- readLines("merged.txt")
mdata <- readLines("merged2.txt")

s1 <- unlist(lapply(strsplit(mdata, " in chromosomes "), "[",2))
l1 <- unlist(lapply(strsplit(s1, " and "), "[", 1))
l2 <- unlist(lapply(strsplit(s1, " and "), "[", 2))

chr1 <- unlist(lapply(strsplit(l1,":"), "[", 1))
posi1 <- unlist(lapply(strsplit(l1,":"), "[", 2))
pos1 <- (as.numeric(unlist(lapply(strsplit(posi1,"-"), "[", 1))) + as.numeric(unlist(lapply(strsplit(posi1,"-"), "[", 2)))) / 2
chr2 <- unlist(lapply(strsplit(l2,":"), "[", 1))
posi2 <- unlist(lapply(strsplit(l2,":"), "[", 2))
pos2 <- (as.numeric(unlist(lapply(strsplit(posi2,"-"), "[", 1))) + as.numeric(unlist(lapply(strsplit(posi2,"-"), "[", 2)))) / 2

plot(pos1[which(chr1 == "1" & chr2 == "1")], pos2[which(chr1 == "1" & chr2 == "1")], pch=20, col=rgb(1,0,0,0.1))
plot(pos1[which(chr1 == "2" & chr2 == "2")], pos2[which(chr1 == "2" & chr2 == "2")], pch=20, col=rgb(1,0,0,0.1))
plot(pos1[which(chr1 == "3" & chr2 == "3")], pos2[which(chr1 == "3" & chr2 == "3")], pch=20, col=rgb(1,0,0,0.1))
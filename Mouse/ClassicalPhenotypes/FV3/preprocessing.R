# Data Pre-processing
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written March, 2009
#

setwd("E:/Mouse/ClassicalPhenotypes/FV3")

crossFF <- read.table("cross_F2_FF_input.csv", header=TRUE, sep=",", colClasses="character")
crossFF[which(crossFF[,"Sex"]==0),"Sex"] <- 'f'
crossFF[which(crossFF[,"Sex"]==1),"Sex"] <- 'm'
crossFF[3:nrow(crossFF), grep("Mit", colnames(crossFF))] <- apply(crossFF[3:nrow(crossFF), grep("Mit", colnames(crossFF))],2,function(x){
  r <- rep(NA, length(x))
  r[which(x==1)] <- "A"
  r[which(x==2)] <- "H"
  r[which(x==3)] <- "B"
  r
})
write.table(crossFF,"cross_F2_FF.csv", quote=FALSE, row.names=FALSE, sep=",")

crossNF <- read.table(file="cross_F2_NF_input.csv", header=TRUE, sep=",", colClasses="character")
crossNF[3:nrow(crossNF), 63:ncol(crossNF)] <- apply(crossNF[3:nrow(crossNF), 63:ncol(crossNF)],2,function(x){
  r <- rep(NA, length(x))
  r[which(x=="1")] <- "A"
  r[which(x=="3")] <- "H"
  r[which(x=="2")] <- "B"
  r
})
write.table(crossNF,"cross_F2_NF.csv", row.names=FALSE, quote=FALSE, sep=",")

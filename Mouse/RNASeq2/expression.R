
# Loading data
setwd("E:/Mouse/RNA/ArrayDesign/Atlas data")
arrays <- read.table("Annotation/arrays.txt", header=TRUE, sep="\t", colClasses="character")
alldata <- read.table("Analysis/geneexpression.txt", sep="\t", header=TRUE)

ofInterest <- c("Lep", "Lepr", "Tph1", "Tph2")

msubset <- alldata[which(alldata[,"mgi_symbol"] %in% ofInterest),]

HT_BFMI <- arrays[which(arrays[,"Tissue"] == "HT" & arrays[,"Strain"] == "BFMI"),"AtlasID"]
HT_B6N  <- arrays[which(arrays[,"Tissue"] == "HT" & arrays[,"Strain"] == "B6N") ,"AtlasID"]

GF_BFMI <- arrays[which(arrays[,"Tissue"] == "GF" & arrays[,"Strain"] == "BFMI"),"AtlasID"]
GF_B6N  <- arrays[which(arrays[,"Tissue"] == "GF" & arrays[,"Strain"] == "B6N") ,"AtlasID"]

msubset <- cbind(msubset, HT_BFMIvsB6N = NA)
msubset <- cbind(msubset, GF_BFMIvsB6N = NA)

for(x in 1:nrow(msubset)){
  msubset[x, "HT_BFMIvsB6N"] <- t.test(msubset[x, HT_BFMI], msubset[x, HT_B6N])$p.value
  msubset[x, "GF_BFMIvsB6N"] <- t.test(msubset[x, GF_BFMI], msubset[x, GF_B6N])$p.value
}

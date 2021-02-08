setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
mdata <- read.csv("filteredSNPs.txt", sep = "\t")

mdata <- mdata[which(mdata[, "CHROM"] == "ENA|CM004567|CM004567.1"),]

ss <- min(mdata[, "POS"])
se <- max(mdata[, "POS"])

ns <- as.numeric(ncol(mdata) - apply(apply(mdata[,10:ncol(mdata)],1,is.na),2,sum))
depth <- as.numeric(gsub("DP=", "", unlist(lapply(strsplit(mdata[, "INFO"], ";"), "[", 1))))
ii <- which(depth / ns >= 10)

mdata <- mdata[ii,]
ns <- as.numeric(ncol(mdata) - apply(apply(mdata[,10:ncol(mdata)],1,is.na),2,sum))
depth <- as.numeric(gsub("DP=", "", unlist(lapply(strsplit(mdata[, "INFO"], ";"), "[", 1))))


ym <- max(depth/ ns)

nssnps <- c(85981710,85981711,85982615,85984154,85987197,85988705,86015278,86015259,86013169,86008103,86008016,86079098,86081790,86081887,86085160,86085714,86089407,86208927,86208939,86208960,86209097,86209263)
syn <- c(85979976,85981703,85982631,85988712,85991559,85993377,85993386,86015270,86206785,86208859,86208883,86208889,86208958)

plot(c(ss, se), c(0,ym), ylab = "Depth", xlab = "Position (bp)", t = 'n')

nsii <- which(mdata[, "POS"] %in% nssnps)
sii <- which(mdata[, "POS"] %in% syn)

points(mdata[, "POS"], depth / ns, t = 'h', col = c("gray"))
points(mdata[nsii, "POS"], (depth / ns)[nsii], t = 'h', col = c("red"), lwd=3)
points(mdata[sii, "POS"], (depth / ns)[sii], t = 'h', col = c("blue"), lwd=2)



length(which(mdata[, "POS"] %in% nssnps))
length(nssnps)
length(which(mdata[, "POS"] %in% syn))
length(syn)

depth[pp] / ns[pp]

setwd("d:/")
mm <- read.table("positions.txt")
setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
mdata <- read.csv("filteredSNPs.txt", sep = "\t")
mdata <- mdata[which(mdata[, "CHROM"] == "ENA|CM004567|CM004567.1"),]

for(x in 1: nrow(mm)){
  i <- which(mdata[, "POS"] == mm[x,1])
  cat(mm[x,1], "\t", mdata[i, "QUAL"], "\n", file = "posQual.txt", append = TRUE)
}
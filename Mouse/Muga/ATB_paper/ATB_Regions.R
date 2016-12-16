# Load the Allele Transmission Biased regions

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
sPat28 <- read.table("Analysis/TransmissionBias_Pat_0.01_28.txt", sep="\t")
sMat28 <- read.table("Analysis/TransmissionBias_Mat_0.01_28.txt", sep="\t")

setwd("D:/Edrive/Mouse/ClassicalPhenotypes/AIL/Analysis")

qtls <- read.table("qtls_d63.txt",sep="\t")

op <- par(mfrow=c(2,1))
plot(qtls[rownames(sPat28),"marker"],main="LOD scores d63, pat TRD regions")
abline(h=threshold <- -log10(0.05 / length(rownames(sPat28))))
plot(qtls[rownames(sMat28),"marker"],main="LOD scores d63, mat TRD regions")
abline(h=threshold <- -log10(0.05 / length(rownames(sMat28))))




# Load the genotype call data
setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
locusxdnaheader <- unlist(strsplit(readLines("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", n=16)[16],","))
locusxdnaR <- read.csv("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", header=FALSE, skip = 22)
colnames(locusxdnaR) <- c("Label","plateWell","Date","oligoPoolId","bundleId", "status", "Type", "Nas", locusxdnaheader[4:length(locusxdnaheader)])
locusxdna <- t(locusxdnaR[which(locusxdnaR[,"Type"] == "calls"),])
locusxdnascores <- t(locusxdnaR[which(locusxdnaR[,"Type"] == "Score_Call"),])
colnames(locusxdna) <- locusxdna[1,]
colnames(locusxdnascores) <- locusxdnascores[1,]
locusxdna <- locusxdna[-c(1:8), ] # Remove some unneeded headers
locusxdnascores <- locusxdnascores[-c(1:8), ] # Remove some unneeded headers
locusxdnascores[which(locusxdna == "U")] <- NA # Use the correct NA values
locusxdna[which(locusxdna == "U")] <- NA # Use the correct NA values

write.table(locusxdnascores, "GCscores.GTS.txt", sep="\t", quote=FALSE)
write.table(locusxdna, "GTS.txt", sep="\t", quote=FALSE)

# Start here
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
setwd("D:/Edrive/Mouse/DNA/MegaMuga/")

genotypes <- read.table("GTS.txt",sep="\t", header=TRUE, na.strings=c(0,"-","NA"), colClasses="character", check.names=FALSE)
scores <- read.table("GCscores.GTS.txt",sep="\t", row.names=1, header=TRUE, na.strings=c(0,"-","NA"), check.names=FALSE)
phenotypes <- read.table("Phenotypes/allPhenotypes.txt",sep="\t", row.names=1, header=TRUE, na.strings=c(0,"-","NA"), check.names=FALSE)

phenotypes <- phenotypes[which(phenotypes[,"Gen."] == 28),]
genotypes <- genotypes[,rownames(phenotypes)]
scores <- scores[,rownames(phenotypes)]

mPsnp <- apply(apply(genotypes,1,is.na),2,sum)
mPind <- apply(apply(genotypes,2,is.na),2,sum)

indTR <- names(which(mPind > 10000))
snpTR <- names(which(mPsnp > 35))

phenotypes <- phenotypes[-which(rownames(phenotypes) %in% indTR),]
genotypes <- genotypes[-which(rownames(genotypes) %in% snpTR),rownames(phenotypes)]
scores <- scores[-which(rownames(scores) %in% snpTR),rownames(phenotypes)]

noSeg <- which(lapply(apply(genotypes,1,table), length) == 1)
genotypes <- genotypes[-noSeg,]
scores <- scores[-noSeg,]

smallG <- which(lapply(apply(genotypes,1,table), min) < 20)
genotypes <- genotypes[-smallG,]
scores <- scores[-smallG,]

lods <- c()
for(x in 1:nrow(genotypes)){
  rN <- anova(lm(phenotypes[, "d70"] ~ phenotypes[, "WG2"] + phenotypes[, "W.Label"] + as.character(genotypes[x,])))
  rW <- anova(lm(phenotypes[, "d70"] ~ phenotypes[, "WG2"] + phenotypes[, "W.Label"] + as.character(genotypes[x,]), weights=as.numeric(scores[x,])))
  lods <- rbind(lods, -log10(c(rN["as.character(genotypes[x, ])", "Pr(>F)"],rW["as.character(genotypes[x, ])", "Pr(>F)"])))
}

plot(lods, xlab="normal", ylab="weighted")
cor(lods)

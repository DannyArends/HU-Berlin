# QTL analysis of the QTL data from Mustapha
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

# TODO use forward selection
library(lme4)
setwd("D:/Edrive/Chicken/ClassicalPhenotypes/QTLanalysis_Mostafa_2010")

genotypes <- read.table("genotypes.txt", sep="\t", na.strings = c("", "."), header=TRUE,colClasses="character")
bodycomposition <- read.table("bodycomposition.txt", sep="\t", na.strings = ".", header=TRUE)[, -1]
growth <- read.table("growth.txt", sep="\t", na.strings = ".", header=TRUE)
map <- read.table("map.txt", sep="\t", row.names=1)

chromosomes <- c(1:29)
names(chromosomes) <- c(1:28,"Z")

chromosomes <- chromosomes[which(names(chromosomes) %in% map[,1])]

map <- map[sort(unlist(lapply(map[,1], function(x){which(names(chromosomes) == x)})), index.return=T)$ix,]

# Duplicated measurements ?!
growth <- growth[-which(duplicated(growth[,"ID"])),]

bodycomposition <- bodycomposition[which(bodycomposition[,"ID"] %in% growth[,"ID"]),]
growth <- growth[which(growth[,"ID"] %in% bodycomposition[,"ID"]),]

rownames(genotypes) <- paste0("C", genotypes[,"ID"])
rownames(growth) <- paste0("C", growth[,"ID"])
rownames(bodycomposition) <- paste0("C", bodycomposition[,"ID"])

# Merge
phenotypes <- cbind(bodycomposition, growth[rownames(bodycomposition),])
phenotypes <- phenotypes[which(rownames(phenotypes) %in% rownames(genotypes)),]
genotypes <- genotypes[rownames(phenotypes),]

markers <- rownames(map)
cofactors <- colnames(phenotypes)[1:7]
phenonames <- colnames(phenotypes)[which(!colnames(phenotypes) %in% cofactors)]
timepoints <- c("BW0", "BW5", "BW10", "BW15", "BW20")

hatch <- as.factor(as.character(phenotypes[,"HatchDate"]))
sire <- as.factor(as.character(phenotypes[,"Sire"]))

for(mname in markers){
  mtable <- table(genotypes[,mname])
  small <- which(mtable < 10)
  genotypes[which(genotypes[,mname] %in% names(small)),mname] <- NA
}

p.LM <- matrix(NA, length(phenonames), length(markers), dimnames=list(phenonames, markers))
p.LMM <- matrix(NA, length(phenonames), length(markers), dimnames=list(phenonames, markers))
for(pname in phenonames){
  for(mname in markers){
    P <- phenotypes[,pname]
    M <- genotypes[,mname]
    mymodel <- lm(P ~  hatch + sire + M)
    p.LM[pname, mname] <- anova(mymodel)["Pr(>F)"]["M",]
    
    hatchLMM <- hatch
    sireLMM <- sire
    naM <- which(is.na(M))
    if(length(naM) > 0){
      P <- P[-naM]
      hatchLMM <- hatchLMM[-naM]
      sireLMM <- sireLMM[-naM]
      M <- M[-naM]
    }
    
    lmm.full <- lmer(P ~  hatchLMM + (1|sireLMM) + M, REML=FALSE)
    lmm.null <- lmer(P ~  hatchLMM + (1|sireLMM), REML=FALSE)
    p.LMM[pname, mname] <- anova(lmm.full, lmm.null)["lmm.full", "Pr(>Chisq)"]
  }
}
lod.LM <- t(-log10(p.LM))
lod.LMM <- t(-log10(p.LMM))
lod.D <- abs(lod.LM - lod.LMM)

plot.data <- function(lods){
  colz <- c("white", colorRampPalette(c("lightblue", "darkblue"))(7))
  breaks <- seq(0,(max(lods)+1), (max(lods)+1) / length(colz))
  op <- par(mar = c(4,15,1,1))
  image(1:nrow(lods),1:ncol(lods), lods, yaxt='n', ylab='', xlab="Marker", col=colz, breaks = breaks)
  axis(2, at=1:ncol(lods), phenonames, las=2)
  box()
  par(xpd=TRUE)
  legend(-30, ncol(lods), apply(round(cbind(breaks, breaks[-1])[-length(breaks),],1),1, paste0, collapse=" - "), fill=colz)
  par(xpd=FALSE)
}

plot.data(lod.LM)
plot.data(lod.LMM)
plot.data(lod.D)


tp.phe <- as.numeric(unlist(phenotypes[,timepoints]))
tp.day <- unlist(lapply(timepoints, rep, nrow(phenotypes)))
tp.date <- rep(phenotypes[,"HatchDate"], length(timepoints))
tp.sire <- rep(phenotypes[,"Sire"], length(timepoints))
tp.ind <- rep(rownames(phenotypes), length(timepoints))

noPhe <- which(is.na(tp.phe))

tp.phe <- tp.phe[-noPhe]
tp.day <- tp.day[-noPhe]
tp.date <- tp.date[-noPhe]
tp.sire <- tp.sire[-noPhe]
tp.ind <- tp.ind[-noPhe]

tp.dayN <- as.numeric(gsub("BW", "", as.character(tp.day)))

intres <- matrix(NA, length(markers), 3, dimnames=list(markers, c("G0","G1", "GxT")))
for(mname in markers){
  M <- genotypes[tp.ind, mname]
  tp.pheLMM <- tp.phe
  tp.dayLMM <- tp.dayN
  tp.dateLMM <- tp.date
  tp.sireLMM <- tp.sire
  tp.indLMM <- tp.ind
  noM <- which(is.na(M))
  if(length(noM) > 0){
    tp.pheLMM <- tp.pheLMM[-noM]
    tp.dayLMM <- tp.dayLMM[-noM]
    tp.dateLMM <- tp.dateLMM[-noM]
    tp.sireLMM <- tp.sireLMM[-noM]
    tp.indLMM <- tp.indLMM[-noM]
    M <- M[-noM]
  }
  
  model.full   <- lmer(tp.pheLMM ~ tp.dateLMM + (1|tp.sireLMM) + tp.dayLMM + (1|tp.indLMM) + M + M:tp.dayLMM, REML=FALSE)
  model.marker <- lmer(tp.pheLMM ~ tp.dateLMM + (1|tp.sireLMM) + tp.dayLMM + (1|tp.indLMM) + M, REML=FALSE)
  model.int    <- lmer(tp.pheLMM ~ tp.dateLMM + (1|tp.sireLMM) + tp.dayLMM + (1|tp.indLMM) + M:tp.dayLMM, REML=FALSE)
  model.null   <- lmer(tp.pheLMM ~ tp.dateLMM + (1|tp.sireLMM) + tp.dayLMM + (1|tp.indLMM), REML=FALSE)
  intres[mname, "G0"] <-  anova(model.marker, model.null)["model.marker", "Pr(>Chisq)"]
  intres[mname, "G1"] <-  anova(model.full, model.int)["model.full", "Pr(>Chisq)"]
  intres[mname, "GxT"] <- anova(model.full, model.marker)["model.full", "Pr(>Chisq)"]
}
lod.Int <- -log10(intres)

op <- par(mfrow=c(2,1))
op <- par(mar = c(4,5,3,1))
plot(c(0, nrow(intres)+1), c(0, max(lod.Int)), t = 'n',xlab="",xaxt='n', ylab="LOD", main="Combined LMM", xaxs="i")
axis(1, at = 1:nrow(intres), markers, las=2, cex.axis=0.8)
s <- 0
for(chr in names(chromosomes)){
  onCHR <- rownames(map[map[,1] == chr,])
  type <- 'l'
  if(length(onCHR) == 1) type = 'p'
  points(s + (1:length(onCHR)), lod.Int[onCHR,"G0"], t = type, lwd=2, pch=18)
  points(s + (1:length(onCHR)), lod.Int[onCHR,"G1"], t = type, lwd=2, col="darkgreen", pch=18)
  points(s + (1:length(onCHR)), lod.Int[onCHR,"GxT"], t = type, col='blue', pch=18)
  s = s + length(onCHR)
  abline(v = s + 0.5)
}
abline(h = -log10(0.05 / nrow(lod.Int)), col = "orange")
abline(h = -log10(0.01 / nrow(lod.Int)), col = "green")

plot(c(0, nrow(intres)+1), c(0, max(lod.Int)), t = 'n',xlab="",xaxt='n', ylab="LOD", main="Timepoint LMM", xaxs="i")
axis(1, at = 1:nrow(intres), markers, las=2, cex.axis=0.8)
s <- 0
for(chr in names(chromosomes)){
  onCHR <- rownames(map[map[,1] == chr,])
  type <- 'l'
  if(length(onCHR) == 1) type = 'p'
  points(s + (1:length(onCHR)),lod.LMM[onCHR,"BW0"], t = type, col='gray', lty=1, pch=18)
  points(s + (1:length(onCHR)),lod.LMM[onCHR,"BW5"], t = type, col='gray', lty=2, pch=18)
  points(s + (1:length(onCHR)),lod.LMM[onCHR,"BW10"], t = type, col='gray', lty=3, pch=18)
  points(s + (1:length(onCHR)),lod.LMM[onCHR,"BW15"], t = type, col='gray', lty=4, pch=18)
  points(s + (1:length(onCHR)),lod.LMM[onCHR,"BW20"], t = type, col='gray', lty=5, pch=18)
  s = s + length(onCHR)
  abline(v = s + 0.5)
}
abline(h = -log10(0.05 / nrow(lod.Int)), col = "orange")
abline(h = -log10(0.01 / nrow(lod.Int)), col = "green")

mname <- "UMA4.034"
M <- genotypes[tp.ind, mname]
tp.pheLMM <- tp.phe
tp.dayLMM <- tp.dayN
tp.dateLMM <- tp.date
tp.sireLMM <- tp.sire
tp.indLMM <- tp.ind
noM <- which(is.na(M))
if(length(noM) > 0){
  tp.pheLMM <- tp.pheLMM[-noM]
  tp.dayLMM <- tp.dayLMM[-noM]
  tp.dateLMM <- tp.dateLMM[-noM]
  tp.sireLMM <- tp.sireLMM[-noM]
  tp.indLMM <- tp.indLMM[-noM]
  M <- M[-noM]
}

model.full   <- lmer(tp.pheLMM ~ tp.dateLMM + tp.sireLMM + tp.dayLMM + (1|tp.indLMM) + M + M:tp.dayLMM, REML=FALSE)
model.marker <- lmer(tp.pheLMM ~ tp.dateLMM + tp.sireLMM + tp.dayLMM + (1|tp.indLMM) + M, REML=FALSE)
anova(model.full, model.marker)

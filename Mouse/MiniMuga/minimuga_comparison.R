setwd("D:/Edrive/Mouse/DNA/MiniMuga")

mdata <- read.csv("miniMUGA-Genotypes.csv", sep=',', na.strings="N", colClasses="character", row.names=1)
mdata <- mdata[mdata[,"Chromosome"] != 0,]
colnames(mdata)[2] <- "Position"
colnames(mdata) <- gsub("BFMI.Gudrun.Brockmann_", "", colnames(mdata))
samples <- rbind(c("M11113", "BFMI860"), c("M28334", "BFMI860-12"), c("M7194", "BFMI861-S1"), c("M7236","BFMI861-S2"))
colnames(mdata)[3:6] <- samples[,2]

annot <- read.csv("miniMUGA-MarkerAnnotations.csv", sep=',', na.strings="N", colClasses="character", row.names=1)
annot <- annot[rownames(mdata),]

different <- which(annot[rownames(mdata),"reference"] != mdata[,"BFMI861-S1"])
difSNPs <- annot[different, c("Chromosome", "Position..b38.")]
allSNPs <- annot[, c("Chromosome", "Position..b38.")]
nrow(difSNPs)
table(annot[which(annot[rownames(mdata),"reference"] != mdata[,"BFMI861-S1"]),"Chromosome"])

mt <- mdata[which(mdata[,"Chromosome"] == "MT"),]
mtgt <- cbind(annot[rownames(mt),], mt)[,c("Chromosome", "Position..b38.", "reference", "alternate", "BFMI860", "BFMI860-12", "BFMI861-S1", "BFMI861-S2")]

mtdiffs <- mtgt[which(lapply(apply(mtgt[, c(3,5:7)],1,table), length) > 1), ]

mtVEP <- mtdiffs[, c("Chromosome", "Position..b38.", "Position..b38.", "reference", "alternate")]
for(x in 1:nrow(mtVEP)){
  cat(unlist(mtVEP[x, ]), 1, "\n")
}
op <- par(mfrow=c(1,2))

mY <- max(as.numeric(difSNPs[,"Position..b38."]))
plot(c(1,20), c(0,mY), t = 'n',xlab="", ylab="Position", xaxt='n', yaxt='n', main="MiniMuga (ALL)")
pos <- 1
for(chr in c(1:19,"X")){
  onChrA <- allSNPs[which(allSNPs[,"Chromosome"] == chr),]
  rect(rep(pos-0.15,nrow(onChrA)) , as.numeric(onChrA[,"Position..b38."]), rep(pos+0.15,nrow(onChrA)), as.numeric(onChrA[,"Position..b38."]))
  rect(pos-0.2, 0, pos+0.2, max(as.numeric(onChrA[,"Position..b38."])))
  pos <- pos + 1
}
axis(1, at = 1:20, paste0("Chr ", c(1:19,"X")), cex=0.8, las=2)
axis(2, at = seq(0, mY, 10000000), seq(0, mY, 10000000) / 1000000, las=2)

mY <- max(as.numeric(difSNPs[,"Position..b38."]))
plot(c(1,20), c(0,mY), t = 'n',xlab="", ylab="Position", xaxt='n', yaxt='n', main="MiniMuga (Informative)")
pos <- 1
for(chr in c(1:19,"X")){
  onChr <- difSNPs[which(difSNPs[,"Chromosome"] == chr),]
  rect(rep(pos-0.15,nrow(onChr)) , as.numeric(onChr[,"Position..b38."]), rep(pos+0.15,nrow(onChr)), as.numeric(onChr[,"Position..b38."]))
  rect(pos-0.2, 0, pos+0.2, max(as.numeric(onChr[,"Position..b38."])))
  pos <- pos + 1
}
axis(1, at = 1:20, paste0("Chr ", c(1:19,"X")), cex=0.8, las=2)
axis(2, at = seq(0, mY, 10000000), seq(0, mY, 10000000) / 1000000, las=2)

snps <- cbind(table(annot[, "Chromosome"])[c(1:19, "X", "Y", "MT")],
      table(annot[which(annot[rownames(mdata),"reference"] != mdata[,"BFMI861-S1"]),"Chromosome"])[c(1:19, "X", "Y", "MT")])

colnames(snps) <- c("All", "Informative")
write.table(snps, file="informativesnps.txt",sep='\t', quote = FALSE)

     
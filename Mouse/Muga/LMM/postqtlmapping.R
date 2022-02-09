# post QTL mapping analysis of the MegaMuga data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
library(parallel)

setwd("D:/Edrive/Mouse/DNA/MegaMuga/inputF2")
genotypes <- read.table(file="cleaned_recoded_genotypes_F2.txt", sep = '\t', check.names = FALSE)
map <- read.table(file="cleaned_map_25MbGap.txt", sep = '\t')
phenotypes <- read.table(file="cleaned_phenotypes_F2.txt", sep = '\t')

timepoints <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")

chrs <- 1:20; names(chrs) <- c(1:19, "X")

chrs.starts <- c(0)
chrs.lengths <- c()
chrs.summed <- 0
chr.gap <- 25000000

for(chr in names(chrs)){
  onChr <- rownames(map[map[, "Chr"] == chr,])
  chr.length <- max(as.numeric(map[onChr, "Mb_NCBI38"]))
  chrs.summed <- chrs.summed + chr.length + chr.gap
  chrs.lengths <- c(chrs.lengths, chr.length + chr.gap)
  chrs.starts <- c(chrs.starts, chrs.summed)
}
chrs.lengths <- c(chrs.lengths, NA)

map[,"Chr"] <- factor(map[,"Chr"], levels = names(chrs))


results.lm <- vector("list", length(timepoints))
names(results.lm) <- timepoints
results.lmm <- vector("list", length(timepoints))
names(results.lmm) <- timepoints

for (tp in timepoints) {
  results.lm[[tp]] <- read.table(file=paste0("QTL/LM_", tp, ".txt"), sep = '\t')
  results.lmm[[tp]] <- read.table(file=paste0("QTL/LMM_", tp, ".txt"), sep = '\t')
}

chr <- 3
onChr <- rownames(map[map[, "Chr"] == chr,])

plot(c(min(map[onChr, "sumPos"]), max(map[onChr, "sumPos"])), c(0, 45), t = 'n', ylab="LOD score", xlab="Position (Mb)", xaxt='n',las=2, cex.axis=1.4, cex.lab=1.4)
points(map[onChr, "sumPos"], -log10(results.lm[["d70"]][onChr,"MF"]), t ='l', col = "black", lty=1)
points(y = rep(-1, length(onChr)), x= map[onChr, "sumPos"], pch ='|', col = "black")

png("QTL/lamdba.lm.png", height=2000, width=700)
op <- par(mfrow=c(length(timepoints),2))
op <- par(cex=0.8)
for (tp in timepoints) {
  expected.FM <- (rank(results.lm[[tp]][,"MF"], ties.method="first")+0.5) / (length(results.lm[[tp]][,"MF"]) + 1)
  maxqtl.FM.obs <- max(-log10(results.lm[[tp]][,"MF"]),na.rm=TRUE)
  maxqtl.FM.exp <- max(-log10(expected.FM),na.rm=TRUE) + 5
  lambdaFM <- round(median(qchisq(1.0 - results.lm[[tp]][,"MF"], 1),na.rm=TRUE) / qchisq(0.5, 1),3)

  plot(c(0,maxqtl.FM.exp), c(0,maxqtl.FM.obs), t = 'n', xlab="Expected LOD", ylab="Observed LOD", main=paste0(tp, " Model, lambda=", lambdaFM))
  points(-log10(expected.FM), -log10(results.lm[[tp]][,"MF"]), xlab="Expected LOD", ylab="Observed LOD", col=as.numeric(map[,"Chr"]), pch=18)
  abline(0, 1)
  
  expected.FMC <- (rank(results.lm[[tp]][,"MF_C"], ties.method="first")+0.5) / (length(results.lm[[tp]][,"MF_C"]) + 1)
  maxqtl.FMC.obs <- max(-log10(results.lm[[tp]][,"MF_C"]),na.rm=TRUE)
  maxqtl.FMC.exp <- max(-log10(expected.FM),na.rm=TRUE) + 5
  lambdaFMC <- round(median(qchisq(1.0 - results.lm[[tp]][,"MF_C"], 1),na.rm=TRUE) / qchisq(0.5, 1),3)
  
  plot(c(0,maxqtl.FMC.exp), c(0,maxqtl.FMC.obs), t = 'n', xlab="Expected LOD", ylab="Observed LOD", main=paste0(tp, " Model (corrected BFMI locus), lambda=", lambdaFMC))
  points(-log10(expected.FMC), -log10(results.lm[[tp]][,"MF_C"]), xlab="Expected LOD", ylab="Observed LOD", col=as.numeric(map[,"Chr"]), pch=18)
  abline(0, 1)
}
dev.off()

png("QTL/lamdba.lmm.png", height=2000, width=700)
op <- par(mfrow=c(length(timepoints),2))
op <- par(cex=0.8)
for (tp in timepoints) {
  expected.FM <- (rank(results.lmm[[tp]][,"Full"], ties.method="first")+0.5) / (length(results.lmm[[tp]][,"Full"]) + 1)
  maxqtl.FM.obs <- max(-log10(results.lmm[[tp]][,"Full"]),na.rm=TRUE)
  maxqtl.FM.exp <- max(-log10(expected.FM),na.rm=TRUE) + 5
  lambdaFM <- round(median(qchisq(1.0 - results.lmm[[tp]][,"Full"], 1),na.rm=TRUE) / qchisq(0.5, 1),3)

  plot(c(0,maxqtl.FM.exp), c(0,maxqtl.FM.obs), t = 'n', xlab="Expected LOD", ylab="Observed LOD", main=paste0(tp, " Model, lambda=", lambdaFM))
  points(-log10(expected.FM), -log10(results.lmm[[tp]][,"Full"]), xlab="Expected LOD", ylab="Observed LOD", col=as.numeric(map[,"Chr"]), pch=18)
  abline(0, 1)
  
  expected.FMC <- (rank(results.lmm[[tp]][,"Full_C"], ties.method="first")+0.5) / (length(results.lmm[[tp]][,"Full_C"]) + 1)
  maxqtl.FMC.obs <- max(-log10(results.lmm[[tp]][,"Full_C"]),na.rm=TRUE)
  maxqtl.FMC.exp <- max(-log10(expected.FM),na.rm=TRUE) + 5
  lambdaFMC <- round(median(qchisq(1.0 - results.lmm[[tp]][,"Full_C"], 1),na.rm=TRUE) / qchisq(0.5, 1),3)
  
  plot(c(0,maxqtl.FMC.exp), c(0,maxqtl.FMC.obs), t = 'n', xlab="Expected LOD", ylab="Observed LOD", main=paste0(tp, " Model (corrected BFMI locus), lambda=", lambdaFMC))
  points(-log10(expected.FMC), -log10(results.lmm[[tp]][,"Full_C"]), xlab="Expected LOD", ylab="Observed LOD", col=as.numeric(map[,"Chr"]), pch=18)
  abline(0, 1)
}
dev.off()
tpcolors <- colorRampPalette(c("lightblue", "darkblue"))(length(timepoints))
tpcolorsC <- colorRampPalette(c("firebrick1", "firebrick4"))(length(timepoints))
names(tpcolors) <- timepoints
names(tpcolorsC) <- timepoints

png("QTL/qtl.lm.png", height=700, width=2000)
plot(c(0, max(map[, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n',las=2, main="LM per timepoint")
for (tp in timepoints) {
  for (chr in names(chrs)) {
    onChr <- rownames(map[map[, "Chr"] == chr,])
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][onChr,"MF"]), t ='l', col = tpcolors[tp])
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][onChr,"MF_C"]), t ='l', col = tpcolorsC[tp])
  }
}
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")
dev.off()


for (tp in timepoints) {
  png(paste0("QTL/qtl.lm_",tp,".png"), height=700, width=2000)
  plot(c(0, max(map[, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n',las=2, main="LM per timepoint")
  for (chr in names(chrs)) {
    onChr <- rownames(map[map[, "Chr"] == chr,])
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][onChr,"MF"]), t ='l', col = tpcolors[tp])
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][onChr,"MF_C"]), t ='l', col = tpcolorsC[tp])
  }
  axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
  abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
  abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")
  dev.off()
}


png("QTL/qtl.lmm.png", height=700, width=2000)
plot(c(0, max(map[, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n',las=2, main="LMM per timepoint")
for (tp in timepoints) {
  for (chr in names(chrs)) {
    onChr <- rownames(map[map[, "Chr"] == chr,])
    points(map[onChr, "sumPos"], -log10(results.lmm[[tp]][onChr,"Full"]), t ='l', col = tpcolors[tp])
    points(map[onChr, "sumPos"], -log10(results.lmm[[tp]][onChr,"Full_C"]), t ='l', col = tpcolorsC[tp])
  }
}
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")
dev.off()

for (tp in timepoints) {
  png(paste0("QTL/qtl.lmm_",tp,".png"), height=700, width=2000)
  plot(c(0, max(map[, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n',las=2, main="LMM per timepoint")
  for (chr in names(chrs)) {
    onChr <- rownames(map[map[, "Chr"] == chr,])
    points(map[onChr, "sumPos"], -log10(results.lmm[[tp]][onChr,"Full"]), t ='l', col = tpcolors[tp])
    points(map[onChr, "sumPos"], -log10(results.lmm[[tp]][onChr,"Full_C"]), t ='l', col = tpcolorsC[tp])
  }
  axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
  abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
  abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")
  dev.off()
}

top.markers <- c("UNC1938399","UNC030576333","JAX00522656","UNC090485124","UNC30294194")

## heatmap of LMM data
matrixLMMpTP <- matrix(unlist(lapply(results.lmm,"[[", 1)),nrow(genotypes), 8)

breaks <- c(0, 2, 5, 9, 15, 100)
heatcols <- c("white", "lightblue", "orange", "darkgreen", "purple")
chrbreaks <- c()
chrsofar <- 0
chrmap <- table(map[,"Chr"])
for(x in names(chrmap)){
  chrsofar <- chrsofar + chrmap[x]
  chrbreaks <- c(chrbreaks, chrsofar)
}

op <- par(mfrow=c(2,1))
op <- par(mar=c(2,12,5,4))
par(xpd=FALSE)
image(1:nrow(genotypes), 1:8, -log10(matrixLMMpTP), breaks=breaks, col = heatcols, xaxt='n', yaxt='n', ylab="Time point", xlab="", main= "Per time point LMM QTL mapping")
axis(2, at=1:8, timepoints, las=2)
axis(1, at=chrbreaks - (.5 * chrmap), names(chrmap), las=1)
box()
abline(v = chrbreaks)
abline(h = 1:8 + 0.5)
onChr <- rownames(map[map[, "Chr"] == 3,])
par(xpd=TRUE)
legend(-2800, 8.5, paste0(breaks[-length(breaks)], " - ",c(breaks[-1])), fill=heatcols, title="LOD score")

pmarker <- which(rownames(map) %in% top.markers)
pmarker[2] <- pmarker[2] - 150
pmarker[3] <- pmarker[3] + 150
arrows(pmarker, 9.5,which(rownames(map) %in% top.markers), 8.5, length=0.1, col=1, lwd=1)
text(pmarker, 10, paste0("nR",1:5),cex=1.0)

## heatmap of LMM corrected data
matrixLMMCpTP <- matrix(unlist(lapply(results.lmm,"[[", 2)),nrow(genotypes), 8)

op <- par(mar=c(2,12,5,4))
par(xpd=FALSE)
image(1:nrow(genotypes), 1:8, -log10(matrixLMMCpTP), breaks=breaks, col = heatcols, xaxt='n', yaxt='n', ylab="Time point", xlab="", main= "Per time point LMM-MQM mapping")
axis(2, at=1:8, timepoints, las=2)
axis(1, at=chrbreaks - (.5 * chrmap), names(chrmap), las=1)
box()
abline(v = chrbreaks)
abline(h = 1:8 + 0.5)
onChr <- rownames(map[map[, "Chr"] == 3,])
par(xpd=TRUE)

pmarker <- which(rownames(map) %in% top.markers)
pmarker[2] <- pmarker[2] - 150
pmarker[3] <- pmarker[3] + 150
arrows(pmarker, 9.5,which(rownames(map) %in% top.markers), 8.5, length=0.1, col=1, lwd=1)
text(pmarker, 10, paste0("nR",1:5),cex=1.0)

#legend(-1500, 8, paste0(breaks[-length(breaks)], " - ",c(breaks[-1])), fill=heatcols, title="LOD score")




global.lmm <- read.table(file="QTL/LMM_Global.txt", sep = '\t')
significant <- rownames(global.lmm[which(-log10(global.lmm[,"F_M_C"]) >= -log10(0.05/1546)),])

cbind(map[significant,], -log10(global.lmm[significant,"F_M_C"]), -log10(global.lmm[significant,"F_N_C"]))

maxqtl <- max(-log10(global.lmm[,"F_M"]),na.rm=TRUE)

chr <- 5
onChr <- rownames(map[map[, "Chr"] == chr,])

plot(c(min(map[onChr, "sumPos"]), max(map[onChr, "sumPos"])), c(0, 10), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n',las=2)
for (chr in "5") {
  onChr <- rownames(map[map[, "Chr"] == chr,])
  #points(map[onChr, "sumPos"], -log10(global.lmm[onChr,"F_M"]), t ='l', col = "black")
  points(map[onChr, "sumPos"], -log10(global.lmm[onChr,"F_M_C"]), t ='l', col = "blue")
  points(map[onChr, "sumPos"], -log10(global.lmm[onChr,"F_M_C"]), t ='p', col = "blue", pch=19, cex=0.5)
}
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/1546), lty=2, col="orange")
abline(h = -log10(0.01/1546), lty=2, col="green")
#legend("topright", c("LMM mapping", "LMM MQM mapping"), col=c("black", "blue"), lwd=1)






par(xpd=FALSE)
chr <- 3
onChr <- rownames(map[map[, "Chr"] == chr,])

jObes1Loc <- 36691603
dropJobes1 <- which(map[onChr,"Mb_NCBI38"] > jObes1Loc - 5000000 & map[onChr,"Mb_NCBI38"] < jObes1Loc + 5000000)

LMMMQM <- -log10(global.lmm[onChr,"F_M_C"])
LMM <- -log10(global.lmm[onChr,"F_M"])
LMMMQM[dropJobes1] <- LMM[dropJobes1]

plot(c(min(map[onChr, "sumPos"]), max(map[onChr, "sumPos"]) - 80000000), c(0, 15), t = 'n', ylab="LOD score", xlab="Position (Mb)", xaxt='n',las=2, cex.axis=1.4, cex.lab=1.4)
points(map[onChr, "sumPos"], LMM, t ='l', col = "black", lty=2)
points(map[onChr, "sumPos"], LMMMQM, t ='l', col = "blue",lwd = 2)
axis(1, at = chrs.starts[chr] + seq(0, chrs.lengths[chr], 10000000), seq(0, chrs.lengths[chr], 10000000) / 1000000, cex.axis=1.4)
abline(h = -log10(0.05/1546), lty=2, col="orange")
abline(h = -log10(0.01/1546), lty=2, col="green")
legend("topright", c("LMM mapping", "LMM MQM mapping"), col=c("black", "blue"), lwd=c(1,2), lty=c(2,1), cex=1.2)
par(xpd=TRUE)
arrows(map["UNC5048297", "sumPos"], 12, map["UNC5048297", "sumPos"], 15, length=0.1, col=1, lwd=1.2)
text(map["UNC5048297", "sumPos"], 11, expression(paste(italic("jObes1"), "QTL")), cex=1.2)
arrows(map["UNC030576333", "sumPos"], 8, map["UNC030576333", "sumPos"], 7, length=0.1, col=1, lwd=1.5)
text(map["UNC030576333", "sumPos"], 8.5, "nR2", cex=1.2)
arrows(map["JAX00522656", "sumPos"], 7, map["JAX00522656", "sumPos"], 6, length=0.1, col=1, lwd=1.5)
text(map["JAX00522656", "sumPos"], 7.5, "nR3", cex=1.2)
par(xpd=FALSE)

above <- rownames(global.lmm)[which(-log10(global.lmm[,"F_M_C"]) >= -log10(0.05/nrow(genotypes)))]
above <- above[which(above %in% rownames(map))]

global.lmm[c("UNC1938399","UNC030576333","JAX00522656","JAX00174587","UNC30294194"),]



signTable <- cbind(map[above,1:2], round(-log10(global.lmm[above, "F_M_C"]),2))

onChr <- rownames(map[map[, "Chr"] == 1,])
cbind(map[onChr,], round(-log10(global.lmm[onChr, "F_M_C"]),2))

onChr <- rownames(map[map[, "Chr"] == 3,])
cbind(map[onChr,], round(-log10(global.lmm[onChr, "F_M_C"]),2))

onChr <- rownames(map[map[, "Chr"] == 5,])
cbind(map[onChr,], round(-log10(global.lmm[onChr, "F_M_C"]),2))

onChr <- rownames(map[map[, "Chr"] == 9,])
cbind(map[onChr,], round(-log10(global.lmm[onChr, "F_M_C"]),2))

onChr <- rownames(map[map[, "Chr"] == 11,])
cbind(map[onChr,], round(-log10(global.lmm[onChr, "F_M_C"]),2))

onChr <- rownames(map[map[, "Chr"] == 19,])
cbind(map[onChr,], round(-log10(global.lmm[onChr, "F_M_C"]),2))

length(which(map[,"Chr"] == chr))
length(rownames(map[map[, "Chr"] == chr,]))

table(unlist(genotypes["UNC1938399",]))
table(unlist(genotypes["UNC030576333",]))
table(unlist(genotypes["JAX00522656",]))
table(unlist(genotypes["UNC9857021",]))
table(unlist(genotypes["UNC10016752",]))
table(unlist(genotypes["JAX00174587",]))
table(unlist(genotypes["UNC19258918",]))
table(unlist(genotypes["UNC30294194",]))



library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

regions <- c("1:149553681:154868088", "3:26989539:35953921", "3:49901885:52973026", 
             "5:103600922:108555646", "5:117544287:120430728", "9:86816288:99363348", 
             "11:23076977:29036920", "19:37825545:40410259")
genesInRegion <- NULL
for(r in regions){
  res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("chromosomal_region", "biotype"), values = list(r, "protein_coding"), mart = bio.mart)
  write.table(res.biomart, file=paste0("QTL/genes_", gsub(":", "-",r), ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
  genesInRegion <- rbind(genesInRegion, res.biomart)
}
genesInRegion <- genesInRegion[which(genesInRegion[,"chromosome_name"] %in% names(chrs)),]
genesInRegion <- genesInRegion[-137,]
rownames(genesInRegion) <- genesInRegion[,"external_gene_name"]

mmu00010 <- read.table("kegg/mmu00010_Glycolysis Gluconeogenesis.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04910 <- read.table("kegg/mmu04910_Insulin signaling pathway.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04940 <- read.table("kegg/mmu04940_Type I diabetes mellitus.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04930 <- read.table("kegg/mmu04930_Type II diabetes mellitus.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04950 <- read.table("kegg/mmu04950_MODY.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04931 <- read.table("kegg/mmu04931_Insulin resistance.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04975 <- read.table("kegg/mmu04975_Fat digestion and absorption.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu00061 <- read.table("kegg/mmu00061_Fatty acid biosynthesis.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu01040 <- read.table("kegg/mmu01040_Biosynthesis of unsaturated fatty acids.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04920 <- read.table("kegg/mmu04920_Adipocytokine signaling pathway.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04923 <- read.table("kegg/mmu04923_Regulation of lipolysis in adipocytes.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04973 <- read.table("kegg/mmu04973_Carbohydrate digestion and absorption.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)
mmu04974 <- read.table("kegg/mmu04974_Protein digestion and absorption.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)

biomart.mmu00010 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu00010), mart = bio.mart)
biomart.mmu04910 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04910), mart = bio.mart)
biomart.mmu04940 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04940), mart = bio.mart)
biomart.mmu04930 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04930), mart = bio.mart)
biomart.mmu04950 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04950), mart = bio.mart)
biomart.mmu04931 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04931), mart = bio.mart)
biomart.mmu04975 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04975), mart = bio.mart)
biomart.mmu00061 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu00061), mart = bio.mart)
biomart.mmu01040 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu01040), mart = bio.mart)
biomart.mmu04920 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04920), mart = bio.mart)
biomart.mmu04923 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04923), mart = bio.mart)
biomart.mmu04973 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04973), mart = bio.mart)
biomart.mmu04974 <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("external_gene_name"), values = rownames(mmu04974), mart = bio.mart)

biomart.mmu00010 <- biomart.mmu00010[which(biomart.mmu00010[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04910 <- biomart.mmu04910[which(biomart.mmu04910[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04940 <- biomart.mmu04940[which(biomart.mmu04940[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04930 <- biomart.mmu04930[which(biomart.mmu04930[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04950 <- biomart.mmu04950[which(biomart.mmu04950[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04931 <- biomart.mmu04931[which(biomart.mmu04931[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04975 <- biomart.mmu04975[which(biomart.mmu04975[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu00061 <- biomart.mmu00061[which(biomart.mmu00061[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu01040 <- biomart.mmu01040[which(biomart.mmu01040[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04920 <- biomart.mmu04920[which(biomart.mmu04920[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04923 <- biomart.mmu04923[which(biomart.mmu04923[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04973 <- biomart.mmu04973[which(biomart.mmu04973[,"chromosome_name"] %in% names(chrs)),]
biomart.mmu04974 <- biomart.mmu04974[which(biomart.mmu04974[,"chromosome_name"] %in% names(chrs)),]


biomart.mmu00010 <- biomart.mmu00010[-61,]

rownames(biomart.mmu00010) <- biomart.mmu00010[,"external_gene_name"]
rownames(biomart.mmu04910) <- biomart.mmu04910[,"external_gene_name"]
rownames(biomart.mmu04940) <- biomart.mmu04940[,"external_gene_name"]
rownames(biomart.mmu04930) <- biomart.mmu04930[,"external_gene_name"]
rownames(biomart.mmu04950) <- biomart.mmu04950[,"external_gene_name"]
rownames(biomart.mmu04931) <- biomart.mmu04931[,"external_gene_name"]
rownames(biomart.mmu04975) <- biomart.mmu04975[,"external_gene_name"]
rownames(biomart.mmu00061) <- biomart.mmu00061[,"external_gene_name"]
rownames(biomart.mmu01040) <- biomart.mmu01040[,"external_gene_name"]
rownames(biomart.mmu04920) <- biomart.mmu04920[,"external_gene_name"]
rownames(biomart.mmu04923) <- biomart.mmu04923[,"external_gene_name"]
rownames(biomart.mmu04973) <- biomart.mmu04973[,"external_gene_name"]

res.mmu00010 <- biomart.mmu00010[which(rownames(biomart.mmu00010) %in% rownames(genesInRegion)),]
res.mmu04910 <- biomart.mmu04910[which(rownames(biomart.mmu04910) %in% rownames(genesInRegion)),]
res.mmu04940 <- biomart.mmu04940[which(rownames(biomart.mmu04940) %in% rownames(genesInRegion)),]
res.mmu04930 <- biomart.mmu04930[which(rownames(biomart.mmu04930) %in% rownames(genesInRegion)),]
res.mmu04950 <- biomart.mmu04950[which(rownames(biomart.mmu04950) %in% rownames(genesInRegion)),]
res.mmu04931 <- biomart.mmu04931[which(rownames(biomart.mmu04931) %in% rownames(genesInRegion)),]
res.mmu04975 <- biomart.mmu04975[which(rownames(biomart.mmu04975) %in% rownames(genesInRegion)),]
res.mmu00061 <- biomart.mmu00061[which(rownames(biomart.mmu00061) %in% rownames(genesInRegion)),]
res.mmu01040 <- biomart.mmu01040[which(rownames(biomart.mmu01040) %in% rownames(genesInRegion)),]
res.mmu04920 <- biomart.mmu04920[which(rownames(biomart.mmu04920) %in% rownames(genesInRegion)),]
res.mmu04923 <- biomart.mmu04923[which(rownames(biomart.mmu04923) %in% rownames(genesInRegion)),]
res.mmu04973 <- biomart.mmu04973[which(rownames(biomart.mmu04973) %in% rownames(genesInRegion)),]

res.mmu00010
res.mmu04910
res.mmu04940
res.mmu04930
res.mmu04950
res.mmu04931
res.mmu04975
res.mmu00061
res.mmu01040
res.mmu04920
res.mmu04923
res.mmu04973


kegg <- cbind(kegg, res.biomart[rownames(kegg),])
colnames(kegg)[1:3] <- c("kegg_id", "kegg_description", "ko_id")

kegg[which(kegg[, "mgi_symbol"] %in% res.all[,"mgi_symbol"]),]



                          
expected.FM <- (rank(global.lmm[,"F_M"], ties.method="first")+0.5) / (length(global.lmm[,"F_M"]) + 1)
expected.FMC <- (rank(global.lmm[,"F_M_C"], ties.method="first")+0.5) / (length(global.lmm[,"F_M_C"]) + 1)

maxqtl.FM <- max(-log10(global.lmm[,"F_M"]),na.rm=TRUE)
maxqtl.FMC <- max(-log10(global.lmm[,"F_M_C"]),na.rm=TRUE)

plot(-log10(expected.FM), -log10(global.lmm[,"F_M"]), xlab="Expected LOD", ylab="Observed LOD")
points(-log10(expected.FMC), -log10(global.lmm[,"F_M_C"]), xlab="Expected LOD", ylab="Observed LOD", col='green')
abline(a=0, b=1)

round(median(qchisq(1.0 - global.lmm[,"F_M"], 1),na.rm=TRUE) / qchisq(0.5, 1),3)
round(median(qchisq(1.0 - global.lmm[,"F_M_C"], 1),na.rm=TRUE) / qchisq(0.5, 1),3)


### plot for region with LOD drop 49,901,885 - 52,973,026
library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
region <- "3:49901885:52973026"

res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("chromosomal_region", "biotype"), values = list(region, "protein_coding"), mart = bio.mart)


nR3region <- rownames(map[map[, "Chr"] == 3 & map[, "Mb_NCBI38"] >= 49901885 & map[, "Mb_NCBI38"] <= 52973026,])
nR3.start <- min(map[nR3region,"Mb_NCBI38"])
nR3.end <- max(map[nR3region,"Mb_NCBI38"])


nBFMI <- unlist(lapply(apply(genotypes[nR3region,],1, table), "[", "B"))
nH <- unlist(lapply(apply(genotypes[nR3region,],1, table), "[", "H"))
nB6N <- unlist(lapply(apply(genotypes[nR3region,],1, table), "[", "N"))
nB6N[is.na(nB6N)] <- 0
nNA <- c()
mTable <- apply(genotypes[nR3region,],1, table, useNA = "always")
for (x in 1: length(mTable)){
  colID <- is.na(names(mTable[[x]]))
  nNA <- c(nNA, mTable[[x]][colID])
}

forPaula <- cbind(position = map[nR3region,"Mb_NCBI38"], nB6N, nH, nBFMI, nNA, lod = -log10(global.lmm[nR3region,c("F_M_C")]))
rownames(forPaula) <- nR3region
write.table(forPaula, "GGplotRegionChr3.txt",sep="\t", quote=FALSE)

maxN <- 344
maxS <- maxN / (max(-log10(global.lmm[nR3region,c("F_M_C")])) + 1)

plot(c(nR3.start, nR3.end), y = c(1, 344), t = 'n', xlab="Position (Mb)", ylab="LOD", yaxt='n', xaxt='n',xaxs='i',yaxs='i', cex.lab=1.4)
polygon(c(forPaula[1,1], forPaula[,1], forPaula[60,1]), c(0, forPaula[,5] + forPaula[,4] + forPaula[,3] + forPaula[,2], 0), col="white", border=NA)
polygon(c(forPaula[1,1], forPaula[,1], forPaula[60,1]), c(0, forPaula[,4] + forPaula[,3] + forPaula[,2], 0), col="#2b8cbe", border=NA)
polygon(c(forPaula[1,1], forPaula[,1], forPaula[60,1]), c(0, forPaula[,3] + forPaula[,2], 0), col="#a6bddb", border=NA)
polygon(c(forPaula[1,1], forPaula[,1], forPaula[60,1]), c(0, forPaula[,3], 0), col="#ece7f2", border=NA)
points(x = map[nR3region,"Mb_NCBI38"], y = maxS * -log10(global.lmm[nR3region,c("F_M_C")]), t ='l', lwd = 3, lty=1, col="blue")
abline(h = maxS * -log10(0.05/1546), col='orange', lty=2, lwd=2)
abline(h = maxS * -log10(0.01/1546), col='green', lty=2, lwd=2)


for(x in 1:nrow(res.biomart)){
  rect(as.numeric(res.biomart[x,"start_position"]), 10 + (25 * x%%4), as.numeric(res.biomart[x,"end_position"]), 20 + (25 * x%%4), border=NA, col="orange")
  text(x = (as.numeric(res.biomart[x,"end_position"]) + as.numeric(res.biomart[x,"start_position"])) / 2,y = 25 + (25 * x%%4), res.biomart[x,"external_gene_name"], cex=1.2)
}


box()
axis(side=1, at = seq(nR3.start, nR3.end, 300000), round(seq(nR3.start, nR3.end, 300000) / 1000000, 1),las=1, cex.axis=1.4)
axis(side=2, at = seq(0, maxN, maxS), seq(0, maxN, maxS)/maxS,las=2, cex.axis=1.4)

legend("bottomright", c("LMM MQM mapping", "5%", "1%", "BFMI", "B6N", "Heterozygous", "Missing"), bg=rgb(1,1,1,0.9),
       border =c(NA,NA,NA,"black","black","black","black"), 
       fill=c(NA,NA,NA, "#2b8cbe", "#a6bddb", "#ece7f2", "white"), 
       col=c("blue", "orange", "green", NA,NA,NA,NA), lwd=c(3,2,2,1,1),lty=c(1,2,2,1,1), pch=c(NA), cex=1.2)

maxN <- 225
maxS <- 225 / (max(-log10(global.lmm[nR3region,c("F_M_C")])) + 1)

par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(c(nR3.start, nR3.end), y = c(1, maxN), t = 'n', xlab="Position (Mb)", ylab="LOD", yaxt='n', xaxt='n')
points(x = map[nR3region,"Mb_NCBI38"], y = maxS * -log10(global.lmm[nR3region,c("F_M_C")]), t ='l', lwd = 3, lty=3)
points(x = map[nR3region,"Mb_NCBI38"], y = nBFMI, t ='l', col="orange")
points(x = map[nR3region,"Mb_NCBI38"], y = nBFMI, t ='p', col="orange", pch=19)
points(x = map[nR3region,"Mb_NCBI38"], y = nH, t ='l', col="gray")
points(x = map[nR3region,"Mb_NCBI38"], y = nH, t ='p', col="gray", pch=19)
points(x = map[nR3region,"Mb_NCBI38"], y = nB6N, t ='l', col="blue")
points(x = map[nR3region,"Mb_NCBI38"], y = nB6N, t ='p', col="blue", pch=19)
points(x = map[nR3region,"Mb_NCBI38"], y = nNA, t ='l', col="red")
points(x = map[nR3region,"Mb_NCBI38"], y = nNA, t ='p', col="red", pch=19)

axis(side=1, at = seq(nR3.start, nR3.end, 300000), round(seq(nR3.start, nR3.end, 300000) / 1000000, 1),las=1)
axis(side=2, at = seq(0, maxN, maxS), seq(0, maxN, maxS)/maxS,las=2)
axis(side=4, at = seq(0, maxN, 50), seq(0, maxN, 50),las=2)
mtext("# Genotypes", side=4, line=3)
legend("topleft", c("QTL", "# BFMI", "# Heterozygous", "# B6N", "# missing"), col=c("black", "orange", "gray", "blue", "red"), lwd=c(3,1,1,1,1),lty=c(3,1,1,1,1), pch=c(NA,19,19,19,19))

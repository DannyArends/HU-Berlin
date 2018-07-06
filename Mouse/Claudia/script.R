#
# Analysis of the eye experiment (Claudia Brockmann)
#
#           Bbs7  Bbs12
#group A:   wt    wt
#group B:   BFMI  wt
#group C:   wt    BFMI
#group D:   BFMI  BFMI


setwd("D:/Edrive/Mouse/RNA/Sequencing/Claudia")

geno <- read.table("Genotypen mRNA.txt", sep = '\t', row.names=1, header=TRUE)

expression.genes <- read.table("Results/expression/results.genes.matrix.FPKM.tsv", sep = "\t", header=TRUE, check.names=FALSE)
expression.iso <- read.table("Results/expression/results.isoforms.matrix.FPKM.tsv", sep = "\t", header=TRUE, check.names=FALSE)

op <- par(mfrow=c(4,2))
op <- par(mar=c(5,10,5,5))
plot(unlist(expression.genes[which(expression.genes[,1] == "ENSMUSG00000037325"),-1])[1:6], col=c(1,1,1,2,2,2,3,3,3,4,4,4), main="Bbs7", ylab="", las=2, pch=18,cex=2, xaxt='n', xlab="")
axis(1, at = c(2, 5), c("Bbs7 (wt)", "Bbs7 (BFMI)"))#, "Bbs7 (wt), Bbs12 (BFMI)", "Bbs7 (BFMI), Bbs12 (BFMI)"))
title(ylab="Expression (RPKM)", mgp=c(5,5,0))
plot(1:10, t = 'n',xlab="", ylab="",xaxt='n', yaxt='n')

isoforms <- rbind(c("Bbs7-201", "ENSMUST00000040148", "(715aa)"),
                  c("Bbs7-202", "ENSMUST00000108155", "(673aa)"),
                  c("Bbs7-203", "ENSMUST00000108156", "(715aa)"),
                  c("Bbs7-204", "ENSMUST00000129671", "(No protein)"),
                  c("Bbs7-205", "ENSMUST00000142333", "(91aa)"),
                  c("Bbs7-206", "ENSMUST00000199136", "(No protein)"))

wtwt <- c("A1-6503-1-A1-1", "A2-6415-2-A2-2", "A3-6501-3-A3-1")
bfmiwt <- c("B1-7385-2-B1-1", "B2-7343-1-B2-1", "B3-7173-2-B2-3")
wtbfmi <- c("C1-6276-1-C1-1", "C2-6291-1-C2-1", "C3-6437-3-C3-2")
bfmibfmi <- c("D1-7321-3-D1-3", "D2-7321-4-D2-1", "D3-7321-7-D3-1")
                  
isoExpr <- NULL
for(x in 1:nrow(isoforms)){
  plot(unlist(expression.iso[which(expression.iso[,1] == isoforms[x,2]),-1])[1:6], col=c(1,1,1,2,2,2,3,3,3,4,4,4), main=paste0(isoforms[x,1], " ", isoforms[x,2], " ", isoforms[x,3]), ylab="", las=2, pch=18, cex=2, xaxt='n', xlab="")
  title(ylab="Expression (RPKM)", mgp=c(5,5,0))
  axis(1, at = c(2, 5), c("Bbs7 (wt)","Bbs7 (BFMI)")) #, , "Bbs7 (wt), Bbs12 (BFMI)", "Bbs7 (BFMI), Bbs12 (BFMI)"))
  isoExpr <- rbind(isoExpr, unlist(expression.iso[which(expression.iso[,1] == isoforms[x,2]),-1]))
}

rownames(isoExpr) <- isoforms[,1]
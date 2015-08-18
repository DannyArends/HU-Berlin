##
# Vulcano plot
##

### We are using: BFMI_RPKM_Qnorm_ANN_AddDom_plusLog2.txt

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/Analysis")

rnaexpression <- read.table("BFMI_RPKM_Qnorm_ANN_AddDom_plusLog2.txt", sep="\t",header=TRUE)
rnaexpression <- rnaexpression[which(rnaexpression[,"Mean.BFMI860.12xB6N.L"] > 0.5 | rnaexpression[,"Mean.B6NxBFMI860.12.L"] > 0.5),]

colz <- as.numeric(abs(log2(rnaexpression[,"Ratio_F1_log"])) > 0.2630344) + as.numeric(rnaexpression[,"tTest_F1_log"] < 0.10)+ as.numeric(rnaexpression[,"tTest_F1_log"] < 0.10) + 1
plot(log2(rnaexpression[,"Ratio_F1_log"]),-log10(rnaexpression[,"tTest_F1_log"]), col=as.numeric(colz ==4) +2,pch=20,ylab="log10(pvalue)",xlab="log2 Fold change")
abline(v=0.2630344)
abline(v=-0.2630344)
abline(h=-log10(0.10))


# t.test < 0.10
#
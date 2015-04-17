#compare

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

RPKMnolog <- read.table("Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep = "\t", header=TRUE, row.names=1)
RPKMlog <- read.table("Analysis/BFMI_RPKM_Qnorm_ANN_AddDomLog2.txt", sep = "\t", header=TRUE, row.names=1)

RPKMcombined <- cbind(RPKMnolog, Mean.BFMI860.12xB6N.L_log = RPKMlog[,"Mean.BFMI860.12xB6N.L"],  Mean.B6NxBFMI860.12.L_log = RPKMlog[,"Mean.B6NxBFMI860.12.L"],  tTest_F1_log = RPKMlog[,"tTest_F1"], Ratio_F1_log = RPKMlog[,"Ratio_F1"], Ratio_F1_Scale_log = RPKMlog[,"Ratio_F1_Scale"],
                                 Mean.BFMI860_log = RPKMlog[,"Mean.BFMI860"], Mean.B6N_log = RPKMlog[,"Mean.B6N"], tTest_Par_log = RPKMlog[,"tTest_Par"], Ratio_Par_log = RPKMlog[,"Ratio_Par"], Ratio_Par_Scale_log = RPKMlog[,"Ratio_Par_Scale"])

write.table(RPKMcombined, file="Analysis/BFMI_RPKM_Qnorm_ANN_AddDom_plusLog2.txt", sep="\t", row.names=FALSE)

# Mean.B6NxBFMI860.12.L = patBFMI
# Mean.BFMI860.12xB6N.L - matBFMI

underT <- RPKMcombined[,"tTest_F1_log"] < 0.1
above1 <- RPKMcombined[,"Mean.B6NxBFMI860.12.L"] > 1 | RPKMcombined[,"Mean.BFMI860.12xB6N.L"] > 1
selected <- which(underT & above1)

S1old <- read.csv("Analysis/S1_DifferentialGeneExpressionF1s.txt", sep = "\t", skip = 11, row.names = 5)

S1log <- RPKMcombined[selected, c("mgi_symbol","chromosome_name", "start_position","mgi_description", 
                                  "Mean.B6NxBFMI860.12.L_log", "Mean.BFMI860.12xB6N.L_log", "Ratio_F1_log", "tTest_F1_log", 
                                  "Mean.B6N_log", "Mean.BFMI860_log", "Ratio_Par_log", "tTest_Par_log", "A.D_B6NxBFMI860.12", "A.D_BFMI860.12xB6N")]

cat("##### LEGEND SUPPLEMENTAL TABLE DIFFERENTIAL GENE EXPRESSION BETWEEN patBFMI AND matBFMI\t(p < 0.1)\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt")
cat("### TOTAL GENE EXPRESSION:\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# patBFMI\t\t\t; BFMI father ,B6N mother\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# matBFMI\t\t\t; B6N father, BFMI mother\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# Mean patBFMI\t\t\t; mean RPKM values (samples were quantile normalized for total reads before calculating RPKM)\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# Ratio_F1\t\t\t; RPKM ratio between F1s, matBFMI / patBFMI\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# Mean B6N\t\t\t; mean RPKM values parental strain B6N\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# Ratio_Par\t\t\t; RPKM ratio between parental strains, BFMI / B6N\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# Add./Dom.\t\t\t; total gene expression compared to parental strains in patBFMI, dominance of one parental strain total gene expression as a trait given as parental strain name, (p < 0.1). If dominance was found, the respective maternal strain's expression in the majority of cases dominates gene expression levels in F1 males\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# TFBS\t\t\t; transcription factor binding sites for sexual hormones estrogen and anrogens, 2000 upstream and 1000 downstream from gene start, specifity/sensitifity >= 90%\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("# \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n", file= "Analysis/S1log.txt", append=TRUE)
cat("mgi_symbol\tChromosome\tGeneStartPosition\tmgi_description\tensembl_gene_id\tMean patBFMI\tMean matBFMI\tRatio_F1\ttTest_F1\tMean B6N\tMean BFMI860\tRatio_Par\ttTest_Par\tAdd./Dom._patBFMI\tAdd./Dom._matBFMI\tTFBS\n", file= "Analysis/S1log.txt", append=TRUE)

S1log[,"Mean.B6NxBFMI860.12.L_log"] <- round(S1log[,"Mean.B6NxBFMI860.12.L_log"], 2)
S1log[,"Mean.BFMI860.12xB6N.L_log"] <- round(S1log[,"Mean.BFMI860.12xB6N.L_log"], 2)
S1log[,"Ratio_F1_log"] <- round(S1log[,"Ratio_F1_log"], 2)
S1log[,"tTest_F1_log"] <- round(S1log[,"tTest_F1_log"], 5)
S1log[,"Mean.B6N_log"] <- round(S1log[,"Mean.B6N_log"], 2)
S1log[,"Mean.BFMI860_log"] <- round(S1log[,"Mean.BFMI860_log"], 2)
S1log[,"Ratio_Par_log"] <- round(S1log[,"Ratio_Par_log"], 2)
S1log[,"tTest_Par_log"] <- round(S1log[,"tTest_Par_log"], 5)

S1log <- S1log[order(S1log[,"tTest_F1_log"]), ]
S1log <- cbind(ensembl_gene_id = rownames(S1log), S1log)

S1log <- cbind(S1log, TFBS = NA)

for(x in 1:nrow(S1log)){
  name <- rownames(S1log)[x]
  idOld <- which(rownames(S1old) == name)
  if(length(idOld) > 0) S1log[x, "TFBS"] <- as.character(S1old[idOld,"TFBS"])
}

#S1log[rownames(S1old), "TFBS"] <- as.character(S1old[,"TFBS"])

S1log <- S1log[, c("mgi_symbol","chromosome_name", "start_position","mgi_description", "ensembl_gene_id",
                                  "Mean.B6NxBFMI860.12.L_log", "Mean.BFMI860.12xB6N.L_log", "Ratio_F1_log", "tTest_F1_log", 
                                  "Mean.B6N_log", "Mean.BFMI860_log", "Ratio_Par_log", "tTest_Par_log", "A.D_B6NxBFMI860.12", "A.D_BFMI860.12xB6N", "TFBS")]

write.table(S1log, file="Analysis/S1log.txt", sep="\t", row.names=FALSE, append=TRUE, col.names=FALSE, quote=FALSE)



S2old <- read.csv("Analysis/S2_ASE_F1s.txt", sep = "\t",skip=40, stringsAsFactors = FALSE)
S2log <- S2old

for(x in 1:nrow(S2log)){
  name <- S2log[x, "ensembl_gene_id"]
  S2log[x, "Mean.matBFMI"] <- round(RPKMlog[name, "Mean.BFMI860.12xB6N.L"],2)
  S2log[x, "Mean.patBFMI"] <- round(RPKMlog[name, "Mean.B6NxBFMI860.12.L"],2)
  S2log[x, "Ratio_F1"] <- round(RPKMlog[name, "Ratio_F1"],2)
  S2log[x, "tTest_F1"] <- round(RPKMlog[name, "tTest_F1"],5)
  S2log[x, "Mean.BFMI860"] <- round(RPKMlog[name, "Mean.BFMI860"],2)
  S2log[x, "Mean.B6N"] <- round(RPKMlog[name, "Mean.B6N"],2)
  S2log[x, "tTest_Par"] <- round(RPKMlog[name, "tTest_Par"],5)
  S2log[x, "Ratio_Par"] <- round(RPKMlog[name, "Ratio_Par"],2)
}

cat(readLines("Analysis/S2_ASE_F1s.txt",n=41), file="Analysis/S2log.txt", sep = "\n")
write.table(S2log, file="Analysis/S2log.txt", sep = "\t", append=TRUE, col.names=FALSE, quote=FALSE, row.names=FALSE)


plot(RPKMnolog[,"tTest_F1"],RPKMlog[,"tTest_F1"])
plot(RPKMnolog[,"Ratio_F1"],RPKMlog[,"Ratio_F1"])

RPKMnolog[which(RPKMnolog[,"mgi_symbol"] == "Peg3"),]
RPKMlog[which(RPKMlog[,"mgi_symbol"] == "Peg3"),]

RPKMnolog[which(RPKMnolog[,"mgi_symbol"] == "H19"),]
RPKMlog[which(RPKMlog[,"mgi_symbol"] == "H19"),]

noLOGBelow <- RPKMnolog[,"tTest_F1"] < 0.005)
LOGBelow <- RPKMlog[,"tTest_F1"] < 0.005)

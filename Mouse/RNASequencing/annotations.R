# annotations.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014
#
# Annotating the RNA sequencing data (pre-processed by MDC)

library(biomaRt)

setwd("E:/Mouse/RNA/Sequencing")

RPKM <- read.csv("2014-07-04_BFMI_RPKM.txt", sep="\t", header=TRUE)                                               # RNA-Seq summary data from MDC
probeannotation <-  RPKM[,1:6]                                                                                    # Split the data from the probe annotation
RPKM <-  RPKM[,-c(1:6)]
colnames(RPKM) <- gsub(".", "-", gsub("RPKM_", "", colnames(RPKM)), fixed=TRUE)                                   # Fix the inconsistencies in the column names

annotation <- read.csv("SampleAnnotation.txt", sep="\t", header=TRUE)                                             # Sample annotation

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")                                                  # Biomart for mouse genes

if(!file.exists("BiomartAnnotation.txt")){
  ENSMUSG <- as.character(probeannotation[,1])                                                                    # Ensemble genes we want to retrieve
  biomartResults <- NULL
  for(x in seq(0, length(ENSMUSG), 1000)){                                                                        # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(ENSMUSG))                                                                       # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")

    res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_id", "mgi_description"), 
                        filters = c("ensembl_gene_id"), 
                        values = ENSMUSG[x:xend], mart = bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
    Sys.sleep(1)
  }
  write.table(biomartResults, file="BiomartAnnotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("BiomartAnnotation.txt", sep="\t", header=TRUE)
}
hasAnnotation <- match(as.character(probeannotation[,1]), biomartResults[,"ensembl_gene_id"])                         # Annotation matching to normdata
biomartResults <- biomartResults[hasAnnotation, ]

AnnotatedData <- cbind(probeannotation, biomartResults)

# Quadriceps
Qcross1 <- c("F1-V-1004_Q", "F1-V-1016_Q", "F1-V-1020_Q")
Qcross2 <- c("F1-V-1000_Q", "F1-V-1008_Q", "F1-V-1012_Q")
QD1 <- RPKM[,Qcross1]                                                                                                 # BFMI cross BFMI860-12xB6N (D1)
QD2 <- RPKM[,Qcross2]                                                                                                 # BFMI cross B6NxBFMI860-12 (D2)
pval <- apply(cbind(QD1,QD2),1,function(x){ return(t.test(x[1:3], x[4:6])$p.value)})                                  # T-test for differences
QF1 <- cbind(QD1, QD2, apply(QD1, 1, mean), apply(QD2, 1, mean), apply(QD1, 1, mean) / apply(QD2, 1, mean), pval)
colnames(QF1)[7:10] <- c("Mean BFMI860-12xB6N", "Mean B6NxBFMI860-12", "Ratio", "tTest")

# Liver
Lcross1 <- c("F1-V-1004_L", "F1-V-1016_L", "F1-V-1020_L")
Lcross2 <- c("F1-V-1000_L", "F1-V-1008_L", "F1-V-1012_L")
LD1 <- RPKM[,Lcross1]                                                                                                 # BFMI cross BFMI860-12xB6N (D1)
LD2 <- RPKM[,Lcross2]                                                                                                 # BFMI cross B6NxBFMI860-12 (D2)
pval <- apply(cbind(LD1,LD2),1,function(x){ return(t.test(x[1:3], x[4:6])$p.value)})                                  # T-test for differences
LF1 <- cbind(LD1, LD2, apply(LD1, 1, mean), apply(LD2, 1, mean), apply(LD1, 1, mean) / apply(LD2, 1, mean), pval)
colnames(LF1)[7:10] <- c("Mean BFMI860-12xB6N", "Mean B6NxBFMI860-12", "Ratio", "tTest")

LBFMIA  <- c("86026522_L", "86026502_L")
LB6NA   <- c("1006954_L", "1006956_L")
LBFMI <- RPKM[,LBFMIA]
LB6N  <- RPKM[,LB6NA]
pval <- apply(cbind(LBFMI,LB6N),1,function(x){ return(t.test(x[1:2], x[3:4])$p.value)})
LP <- cbind(LBFMI, LB6N, apply(LBFMI, 1, mean), apply(LB6N, 1, mean), apply(LBFMI, 1, mean) / apply(LB6N, 1, mean), pval)
colnames(LP)[5:8] <- c("Mean BFMI860", "Mean B6N", "Ratio", "tTest")

updatedData <- cbind(AnnotatedData[,c("g", "mgi_id", "symbol","chr", "start_position", "end_position", "strand", "biotype", "mgi_description")], QF1, LF1, LP)
colnames(updatedData)[1] <- "ensembl_gene_id"

biomartmissing <- which(is.na(updatedData[,"mgi_description"]))
updatedData[,"mgi_description"] <- as.character(updatedData[,"mgi_description"])
updatedData[biomartmissing, "mgi_description"] <- as.character(probeannotation[biomartmissing, "description"])

write.table(updatedData, file="BFMI_RPKM_ANN.txt", sep="\t", row.names=FALSE)

# Quality control of the data

RPKM_MEAN <- apply(RPKM, 1, mean)

notExpressed <- which(RPKM_MEAN <= 1)
RPKM <- RPKM[-notExpressed, ]
cat("Filtered out:", length(notExpressed), "genes (not expressed)\n")

QD1   <- RPKM[,Qcross1]  ; QD2   <- RPKM[,Qcross2]     # Create new subsets, because we reduced the amount of data
LD1   <- RPKM[,Lcross1]  ; LD2   <- RPKM[,Lcross2]
LBFMI <- RPKM[,LBFMIA]   ; LB6N  <- RPKM[,LB6NA]

# Subset the data
findGroups <- function(BFMI, B6N, F1, pExp=0.1, p = 0.05){
  # Find differentially expressed genes between parents
  pval_B6N_E_BFMI <- apply(cbind(BFMI, B6N), 1, function(x){ return(t.test(x[1:2], x[3:4])$p.value)})
  B6N_D_BFMI      <- which(pval_B6N_E_BFMI < pExp)
  cat("Found", length(B6N_D_BFMI), "differentially expressed genes between B6N and BFMI (p <", pExp, ")\n")

  BFMI_B6N    <- cbind(BFMI, B6N)[B6N_D_BFMI,]
  B6N_F1      <- cbind(B6N,  F1)[B6N_D_BFMI,]
  BFMI_F1     <- cbind(BFMI, F1)[B6N_D_BFMI,]
  B6N_F1_BFMI <- cbind(B6N,  F1, BFMI)[B6N_D_BFMI,]

  # Directional T-testing to see what is happening
  pval_B6N_L_BFMI <- apply(BFMI_B6N, 1, function(x){ return(t.test(x[1:2], x[3:4], alternative="less")$p.value)})        # B6N  > BFMI
  cat("Found", sum(pval_B6N_L_BFMI < 0.1), "B6N > BFMI\n")
  pval_B6N_S_BFMI <- apply(BFMI_B6N, 1, function(x){ return(t.test(x[1:2], x[3:4], alternative="greater")$p.value)})     # B6N  < BFMI
  cat("Found", sum(pval_B6N_S_BFMI < 0.1), "B6N < BFMI\n")

  pval_B6N_L  <- apply(B6N_F1 ,  1, function(x){ return(t.test(x[1:2], x[3:5], alternative="greater")$p.value)})         # B6N  > LD1
  pval_B6N_S  <- apply(B6N_F1 ,  1, function(x){ return(t.test(x[1:2], x[3:5], alternative="less")$p.value)})            # B6N  < LD1
  pval_BFMI_L <- apply(BFMI_F1,  1, function(x){ return(t.test(x[1:2], x[3:5], alternative="greater")$p.value)})         # BFMI > LD1
  pval_BFMI_S <- apply(BFMI_F1 , 1, function(x){ return(t.test(x[1:2], x[3:5], alternative="less")$p.value)})            # B6N  < LD1

  # Test LD1 for Additive/Dominance or Undetermined
  G1  <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S <  p & pval_BFMI_L <  p; wG1 <- which(G1);   cat("Group  1: B6N  < D1/D2  < BFMI:", length(wG1), "\n")
  G2  <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S >= p & pval_BFMI_L <  p; wG2 <- which(G2);   cat("Group  2: B6N == D1/D2  < BFMI:", length(wG2), "\n")
  G3  <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S <  p & pval_BFMI_L >= p; wG3 <- which(G3);   cat("Group  3: B6N  < D1/D2 == BFMI:", length(wG3), "\n")
  GU1 <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S >= p & pval_BFMI_L >= p; wGU1 <- which(GU1); cat("Group U1: B6N == D1/D2 == BFMI:", length(which(GU1)), "\n")

  G4  <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L <  p & pval_BFMI_S <  p; wG4 <- which(G4);   cat("Group  4: BFMI  < D1/D2  < B6N:", length(wG4), "\n")
  G5  <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L >= p & pval_BFMI_S <  p; wG5 <- which(G5);   cat("Group  5: BFMI == D1/D2  < B6N:", length(wG5), "\n")
  G6  <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L <  p & pval_BFMI_S >= p; wG6 <- which(G6);   cat("Group  6: BFMI  < D1/D2 == B6N:", length(wG6), "\n")
  GU2 <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L >= p & pval_BFMI_S >= p; wGU2 <- which(GU2); cat("Group U2: BFMI == D1/D2 == B6N:", length(which(GU2)), "\n")
  return(list(B6N_F1_BFMI[G1,], B6N_F1_BFMI[G2,], B6N_F1_BFMI[G3,], B6N_F1_BFMI[GU1,], B6N_F1_BFMI[G4,], B6N_F1_BFMI[G5,], B6N_F1_BFMI[G6,], B6N_F1_BFMI[GU2,],
              list(pval_B6N_S_BFMI, pval_B6N_L_BFMI, pval_B6N_S, pval_B6N_L, pval_BFMI_S, pval_BFMI_L)))
}

updatedData <- cbind(updatedData, "A/D_BFMI860-12xB6N" = rep(NA, nrow(updatedData)))

analysis_LD1 <- findGroups(LBFMI, LB6N, LD1)
updatedData[rownames(LD1), "A/D_BFMI860-12xB6N"] <- "-"
updatedData[rownames(analysis_LD1[[1]]), "A/D_BFMI860-12xB6N"] <- "ADDITIVE"
updatedData[rownames(analysis_LD1[[2]]), "A/D_BFMI860-12xB6N"] <- "B6N"
updatedData[rownames(analysis_LD1[[3]]), "A/D_BFMI860-12xB6N"] <- "BFMI"
updatedData[rownames(analysis_LD1[[4]]), "A/D_BFMI860-12xB6N"] <- "UNKNOWN"
updatedData[rownames(analysis_LD1[[5]]), "A/D_BFMI860-12xB6N"] <- "ADDITIVE"
updatedData[rownames(analysis_LD1[[6]]), "A/D_BFMI860-12xB6N"] <- "BFMI"
updatedData[rownames(analysis_LD1[[7]]), "A/D_BFMI860-12xB6N"] <- "B6N"
updatedData[rownames(analysis_LD1[[8]]), "A/D_BFMI860-12xB6N"] <- "UNKNOWN"

updatedData <- cbind(updatedData, "A/D_B6NxBFMI860-12" = rep(NA, nrow(updatedData)))
analysis_LD2 <- findGroups(LBFMI, LB6N, LD2)
updatedData[rownames(LD2), "A/D_B6NxBFMI860-12"] <- "-"
updatedData[rownames(analysis_LD2[[1]]), "A/D_B6NxBFMI860-12"] <- "ADDITIVE"
updatedData[rownames(analysis_LD2[[2]]), "A/D_B6NxBFMI860-12"] <- "B6N"
updatedData[rownames(analysis_LD2[[3]]), "A/D_B6NxBFMI860-12"] <- "BFMI"
updatedData[rownames(analysis_LD2[[4]]), "A/D_B6NxBFMI860-12"] <- "UNKNOWN"
updatedData[rownames(analysis_LD2[[5]]), "A/D_B6NxBFMI860-12"] <- "ADDITIVE"
updatedData[rownames(analysis_LD2[[6]]), "A/D_B6NxBFMI860-12"] <- "BFMI"
updatedData[rownames(analysis_LD2[[7]]), "A/D_B6NxBFMI860-12"] <- "B6N"
updatedData[rownames(analysis_LD2[[8]]), "A/D_B6NxBFMI860-12"] <- "UNKNOWN"

write.table(updatedData, file="BFMI_RPKM_ANN_AddDom.txt", sep="\t", row.names=FALSE)

# annotations.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014
#
# Annotating the RNA sequencing data done on our own server

library(biomaRt)

setwd("E:/Mouse/DNA/Sequencing/BFMI")
snps <- read.csv("860v2.pes.s.rmd.rg.realigned.recal.haplotclld.Q30DP3HomFilt_NovelSNPer_detailed.out", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors = FALSE) # DNA-Seq SNP data
snps <- snps[-which(snps$FuncClass == "SYNONYMOUS_CODING"), ]                 # Throw away the SYNONYMOUS_CODING functional class

codingRegion  <- unique(snps[grepl("CODING_REGION", snps$Region),"Gene_ID"])
domainRegion  <- unique(snps[!is.na(snps$DomainRegion),"Gene_ID"])
UTR5          <- unique(snps[grepl("5PRIME_UTR", snps$Region),"Gene_ID"])
UTR3          <- unique(snps[grepl("3PRIME_UTR", snps$Region),"Gene_ID"])
SPLICE        <- unique(snps[grepl("SPLICE_SITE", snps$Region),"Gene_ID"])
STARTSTOP     <- unique(snps[grepl("START", snps$FuncClass),"Gene_ID"])
STARTSTOP     <- unique(c(STARTSTOP, unique(snps[grepl("STOP", snps$FuncClass),"Gene_ID"])))

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

RPKM <- read.csv("Analysis/RPKMnorm.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors = FALSE)                                # RNA-Seq RPKM data
annotation <- read.csv("FastQ/sampledescription.txt", sep="\t", header=TRUE)                                                             # Sample annotation
RPKM <- RPKM[, -2]                                                                                                                       # Throw away the second column
colnamesShort <- unlist(lapply(strsplit(colnames(RPKM),"_"),"[",1))

colnames(RPKM) <- annotation[match(colnamesShort, annotation[,"Lib_id"]),"our_name"]                                                     # Correct column names

if(!file.exists("bioMart/BiomartAnnotation.txt")){
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")                                                # Biomart for mouse genes
  ENSMUSG <- as.character(rownames(RPKM))                                                                         # Ensemble genes we want to retrieve
  biomartResults <- NULL
  for(x in seq(0, length(ENSMUSG), 1000)){                                                                        # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(ENSMUSG))                                                                       # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")

    res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("ensembl_gene_id"), 
                        values = ENSMUSG[x:xend], mart = bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
    Sys.sleep(1)
  }
  write.table(biomartResults, file="bioMart/BiomartAnnotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("bioMart/BiomartAnnotation.txt", sep="\t", header=TRUE, colClasses=c("character"))
}

RPKM <- RPKM[,-c(11,12)]                                                                                              # Remove the quadriceps samples
RPKM <- RPKM[-c(28048, 32151, 35810), ]                                                                               # Remove the essentially constant error in t.test

RPKM[is.na(RPKM)] <- 0                                                                                                # Reset the NA's to 0

RPKM <- apply(RPKM, 2, function(x){log(x+1)})


which(rownames(RPKM) == "ENSMUSG00000091771")

# Liver
Lcross1 <- c("F1-V-1004_L", "F1-V-1016_L", "F1-V-1020_L")
Lcross2 <- c("F1-V-1000_L", "F1-V-1008_L", "F1-V-1012_L")
LD1 <- RPKM[,Lcross1]                                                                                                 # BFMI cross BFMI860-12xB6N (D1)
LD2 <- RPKM[,Lcross2]                                                                                                 # BFMI cross B6NxBFMI860-12 (D2)
pval <- apply(cbind(LD1,LD2),1,function(x){ return(t.test(x[1:3], x[4:6])$p.value)})                                  # T-test for differences
LD1mean <- apply(LD1, 1, mean)
LD2mean <- apply(LD2, 1, mean)
ratios  <-  LD1mean / LD2mean
ratios[which(is.infinite(ratios))] <- NA
ratiosScale <- ratios
ratiosToFlip <- which(LD1mean < LD2mean)
ratiosScale[ratiosToFlip] <- -((LD2mean/LD1mean)[ratiosToFlip])
ratiosScale[which(is.infinite(ratiosScale))] <- NA
LF1 <- cbind(LD1, LD2, LD1mean, LD2mean, ratios, ratiosScale, pval)
colnames(LF1)[7:11] <- c("Mean BFMI860-12xB6N L", "Mean B6NxBFMI860-12 L", "Ratio_F1", "Ratio_F1_Scale", "tTest_F1")

LBFMIA  <- c("86026522_L", "86026502_L")
LB6NA   <- c("1006954_L", "1006956_L")
LBFMI <- RPKM[,LBFMIA]
LB6N  <- RPKM[,LB6NA]
pval <- apply(cbind(LBFMI,LB6N),1,function(x){ return(t.test(x[1:2], x[3:4])$p.value)})
LBFMImean <- apply(LBFMI, 1, mean)
LB6Nmean <- apply(LB6N, 1, mean)
ratios  <-  LBFMImean / LB6Nmean
ratios[which(is.infinite(ratios))] <- NA
ratiosScale <- ratios
ratiosToFlip <- which(LBFMImean < LB6Nmean)
ratiosScale[ratiosToFlip] <- -((LB6Nmean/LBFMImean)[ratiosToFlip])
ratiosScale[which(is.infinite(ratiosScale))] <- NA
LP <- cbind(LBFMI, LB6N, LBFMImean, LB6Nmean, ratios, ratiosScale, pval)
colnames(LP)[5:9] <- c("Mean BFMI860", "Mean B6N", "Ratio_Par", "Ratio_Par_Scale", "tTest_Par")

updatedData <- cbind(LF1, LP)

# Subset the data
findGroups <- function(BFMI, B6N, F1, pExp=0.1, p = 0.1){
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
  pval_BFMI_S <- apply(BFMI_F1 , 1, function(x){ return(t.test(x[1:2], x[3:5], alternative="less")$p.value)})            # GFMI < LD1

  # Test LD1 for Additive/Dominance or Undetermined
  G1  <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S <  p & pval_BFMI_L <  p; wG1 <- which(G1);   cat("Group  1: B6N  < D1/D2  < BFMI:", length(wG1), "\n")
  G2  <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S >= p & pval_BFMI_L <  p; wG2 <- which(G2);   cat("Group  2: B6N == D1/D2  < BFMI:", length(wG2), "\n")
  G3  <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S <  p & pval_BFMI_L >= p; wG3 <- which(G3);   cat("Group  3: B6N  < D1/D2 == BFMI:", length(wG3), "\n")
  GU1 <-  pval_B6N_S_BFMI < 0.1 & pval_B6N_S >= p & pval_BFMI_L >= p; wGU1 <- which(GU1); cat("Group U1: B6N == D1/D2 == BFMI:", length(which(GU1)), "\n")

  G4  <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L <  p & pval_BFMI_S <  p; wG4 <- which(G4);   cat("Group  4: BFMI  < D1/D2  < B6N:", length(wG4), "\n")
  G5  <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L <  p & pval_BFMI_S >= p; wG5 <- which(G5);   cat("Group  5: BFMI == D1/D2  < B6N:", length(wG5), "\n")
  G6  <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L >= p & pval_BFMI_S <  p; wG6 <- which(G6);   cat("Group  6: BFMI  < D1/D2 == B6N:", length(wG6), "\n")
  GU2 <-  pval_B6N_L_BFMI < 0.1 & pval_B6N_L >= p & pval_BFMI_S >= p; wGU2 <- which(GU2); cat("Group U2: BFMI == D1/D2 == B6N:", length(which(GU2)), "\n")
  return(list(B6N_F1_BFMI[wG1,], B6N_F1_BFMI[wG2,], B6N_F1_BFMI[wG3,], B6N_F1_BFMI[wGU1,], B6N_F1_BFMI[wG4,], B6N_F1_BFMI[wG5,], B6N_F1_BFMI[wG6,], B6N_F1_BFMI[wGU2,],
              list(pval_B6N_S_BFMI, pval_B6N_L_BFMI, pval_B6N_S, pval_B6N_L, pval_BFMI_S, pval_BFMI_L)))
}

updatedData <- cbind(updatedData, "A/D_BFMI860-12xB6N" = rep(NA, nrow(updatedData)))
updatedData <- cbind(updatedData, "A/D_B6NxBFMI860-12" = rep(NA, nrow(updatedData)))
updatedData <- cbind(updatedData, "E_A_BFMI860-12xB6N_B6N" = rep(NA, nrow(updatedData)))
updatedData <- cbind(updatedData, "E_D_BFMI860-12xB6N_B6N" = rep(NA, nrow(updatedData)))
updatedData <- cbind(updatedData, "E_A_B6NxBFMI860-12_B6N" = rep(NA, nrow(updatedData)))
updatedData <- cbind(updatedData, "E_D_B6NxBFMI860-12_B6N" = rep(NA, nrow(updatedData)))

analysis_LD1 <- findGroups(LBFMI, LB6N, LD1)
updatedData[rownames(LD1), "A/D_BFMI860-12xB6N"] <- "-"
updatedData[(rownames(LD2)), "E_A_BFMI860-12xB6N_B6N"] <- "-"
updatedData[(rownames(LD2)), "E_D_BFMI860-12xB6N_B6N"] <- "-"
alleffectsD1 <- as.numeric(updatedData[,"Mean BFMI860-12xB6N L"]) - as.numeric(updatedData[,"Mean B6N"])
names(alleffectsD1) <- rownames(updatedData)

updatedData[(rownames(analysis_LD1[[1]])), "A/D_BFMI860-12xB6N"] <- "ADDITIVE"
updatedData[(rownames(analysis_LD1[[1]])), "E_A_BFMI860-12xB6N_B6N"] <- alleffectsD1[(rownames(analysis_LD1[[1]]))]
updatedData[(rownames(analysis_LD1[[2]])), "A/D_BFMI860-12xB6N"] <- "B6N"
updatedData[(rownames(analysis_LD1[[2]])), "E_D_BFMI860-12xB6N_B6N"] <- alleffectsD1[(rownames(analysis_LD1[[2]]))]
updatedData[(rownames(analysis_LD1[[3]])), "A/D_BFMI860-12xB6N"] <- "BFMI"
updatedData[(rownames(analysis_LD1[[3]])), "E_D_BFMI860-12xB6N_B6N"] <- alleffectsD1[(rownames(analysis_LD1[[3]]))]
updatedData[(rownames(analysis_LD1[[4]])), "A/D_BFMI860-12xB6N"] <- "UNKNOWN"
updatedData[(rownames(analysis_LD1[[5]])), "A/D_BFMI860-12xB6N"] <- "ADDITIVE"
updatedData[(rownames(analysis_LD1[[5]])), "E_A_BFMI860-12xB6N_B6N"] <- alleffectsD1[(rownames(analysis_LD1[[5]]))]
updatedData[(rownames(analysis_LD1[[6]])), "A/D_BFMI860-12xB6N"] <- "BFMI"
updatedData[(rownames(analysis_LD1[[6]])), "E_D_BFMI860-12xB6N_B6N"] <- alleffectsD1[(rownames(analysis_LD1[[6]]))]
updatedData[(rownames(analysis_LD1[[7]])), "A/D_BFMI860-12xB6N"] <- "B6N"
updatedData[(rownames(analysis_LD1[[7]])), "E_D_BFMI860-12xB6N_B6N"] <- alleffectsD1[(rownames(analysis_LD1[[7]]))]
updatedData[(rownames(analysis_LD1[[8]])), "A/D_BFMI860-12xB6N"] <- "UNKNOWN"

analysis_LD2 <- findGroups(LBFMI, LB6N, LD2)
updatedData[(rownames(LD2)), "A/D_B6NxBFMI860-12"] <- "-"
updatedData[(rownames(LD2)), "E_A_B6NxBFMI860-12_B6N"] <- "-"
updatedData[(rownames(LD2)), "E_D_B6NxBFMI860-12_B6N"] <- "-"
alleffectsD2 <- as.numeric(updatedData[,"Mean B6NxBFMI860-12 L"]) - as.numeric(updatedData[,"Mean B6N"])
names(alleffectsD2) <- rownames(updatedData)

updatedData[(rownames(analysis_LD2[[1]])), "A/D_B6NxBFMI860-12"] <- "ADDITIVE"
updatedData[(rownames(analysis_LD2[[1]])), "E_A_B6NxBFMI860-12_B6N"] <- alleffectsD2[(rownames(analysis_LD2[[1]]))]
updatedData[(rownames(analysis_LD2[[2]])), "A/D_B6NxBFMI860-12"] <- "B6N"
updatedData[(rownames(analysis_LD2[[2]])), "E_D_B6NxBFMI860-12_B6N"] <- alleffectsD2[(rownames(analysis_LD2[[2]]))]
updatedData[(rownames(analysis_LD2[[3]])), "A/D_B6NxBFMI860-12"] <- "BFMI"
updatedData[(rownames(analysis_LD2[[3]])), "E_D_B6NxBFMI860-12_B6N"] <- alleffectsD2[(rownames(analysis_LD2[[3]]))]
updatedData[(rownames(analysis_LD2[[4]])), "A/D_B6NxBFMI860-12"] <- "UNKNOWN"
updatedData[(rownames(analysis_LD2[[5]])), "A/D_B6NxBFMI860-12"] <- "ADDITIVE"
updatedData[(rownames(analysis_LD2[[5]])), "E_A_B6NxBFMI860-12_B6N"] <- alleffectsD2[(rownames(analysis_LD2[[5]]))]
updatedData[(rownames(analysis_LD2[[6]])), "A/D_B6NxBFMI860-12"] <- "BFMI"
updatedData[(rownames(analysis_LD2[[6]])), "E_D_B6NxBFMI860-12_B6N"] <- alleffectsD2[(rownames(analysis_LD2[[6]]))]
updatedData[(rownames(analysis_LD2[[7]])), "A/D_B6NxBFMI860-12"] <- "B6N"
updatedData[(rownames(analysis_LD2[[7]])), "E_D_B6NxBFMI860-12_B6N"] <- alleffectsD2[(rownames(analysis_LD2[[7]]))]
updatedData[(rownames(analysis_LD2[[8]])), "A/D_B6NxBFMI860-12"] <- "UNKNOWN"

noAnnotationData <- which(!rownames(updatedData) %in% biomartResults[,"ensembl_gene_id"])

EmptyMBiomart <- matrix("", length(rownames(updatedData)[noAnnotationData]), ncol(biomartResults))
colnames(EmptyMBiomart) <- colnames(biomartResults)
EmptyMBiomart[,"ensembl_gene_id"] <- rownames(updatedData)[noAnnotationData]

geneAnnotation <- rbind(EmptyMBiomart, biomartResults)

idxes <- match(rownames(updatedData), geneAnnotation[,"ensembl_gene_id"])
geneAnnotation <- geneAnnotation[idxes,]

updatedData <- cbind(geneAnnotation, updatedData)

updatedData <- cbind(updatedData, "CodingRegion"  = rep(0, nrow(updatedData)))
updatedData[which(updatedData$ensembl_gene_id %in% codingRegion), "CodingRegion"] <- 1
updatedData <- cbind(updatedData, "DomainRegion"  = rep(0, nrow(updatedData)))
updatedData[which(updatedData$ensembl_gene_id %in% domainRegion), "DomainRegion"] <- 1
updatedData <- cbind(updatedData, "5' UTR"        = rep(0, nrow(updatedData)))
updatedData[which(updatedData$ensembl_gene_id %in% UTR5), "5' UTR"] <- 1
updatedData <- cbind(updatedData, "3' UTR"        = rep(0, nrow(updatedData)))
updatedData[which(updatedData$ensembl_gene_id %in% UTR3), "3' UTR"] <- 1
updatedData <- cbind(updatedData, "SpliceSite"    = rep(0, nrow(updatedData)))
updatedData[which(updatedData$ensembl_gene_id %in% SPLICE), "SpliceSite"] <- 1
updatedData <- cbind(updatedData, "Stop/Start"    = rep(0, nrow(updatedData)))
updatedData[which(updatedData$ensembl_gene_id %in% STARTSTOP), "Stop/Start"] <- 1

write.table(updatedData, file="Analysis/BFMI_RPKM_Qnorm_ANN_AddDomLog2.txt", sep="\t", row.names=FALSE)


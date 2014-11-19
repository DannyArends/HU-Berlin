# snpToGene.R - Analyze the SNPs and indels called by the GenomeAnalysisToolKit, and summarize them into genes
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Sep, 2014

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
matB6Nsnps  <- read.csv("maternalB6snps_10reads.txt", sep="\t", header=TRUE)                                        # SNPs detected in the maternal B6N F1 cross
matBFMIsnps <- read.csv("maternalBFMIsnps_10reads.txt", sep="\t", header=TRUE)                                      # SNPs detected in the maternal BFMI F1 cross
RPKM        <- read.csv("Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")   # RPKM values from RNA-Seq

GTF <- read.table("GTF/Mus_musculus.GRCm38.76.gtf", sep="\t")                                                       # Gene models
EXONS <- GTF[which(GTF[,3]=="exon"),]

datastr <- strsplit(as.character(EXONS[,9]), "; ")
lengths <- unlist(lapply(datastr, length))

addInfo <- matrix(NA, length(datastr), 2)
for(x in 1:length(datastr)){ addInfo[x, ] <- c(strsplit(datastr[[x]][1]," ")[[1]][2], strsplit(datastr[[x]][lengths[x]]," ")[[1]][2]); }

EXONS <- cbind(addInfo, EXONS)
uniqueGenes <- unique(EXONS[,1])

geneExonsCoupling <- vector("list", length(uniqueGenes))                                                        # Create a list which enumerate the exons per gene
x <- 1
for(gene in uniqueGenes){
  geneExons <- which(EXONS[,1] == gene)
  possibleExons <- EXONS[geneExons[!duplicated(EXONS[geneExons,2])], -11]
  geneExonsCoupling[[x]] <- possibleExons
  cat(x, "Gene", gene, "has", nrow(geneExonsCoupling[[x]]), "/", length(geneExons),"Exons\n")
  x <- x + 1
}

summarizeImprintedSNPsInGene <- function(uniqueGenes, geneExonsCoupling, direction){
  geneSNPCoupling <- vector("list", length(uniqueGenes))
  x <- 1
  for(gene in geneExonsCoupling){                                                                               # Summarize all SNPs from a direction per gene
    onChr <- as.character(direction[,"Chr"]) == as.character(unique(gene[,3]))
    for(exon in 1:nrow(gene)){
      exonStart <- as.numeric(gene[exon,6])
      exonEnd <- as.numeric(gene[exon,7])
      inEXON <- which(direction[onChr,"Loc"] >= exonStart & direction[onChr,"Loc"] <= exonEnd)
      if(length(inEXON) > 0){
        geneSNPCoupling[[x]] <- rbind(geneSNPCoupling[[x]], cbind(direction[onChr,][inEXON,], Exon=exon))
        geneSNPCoupling[[x]] <- geneSNPCoupling[[x]][!duplicated(geneSNPCoupling[[x]][,"ID"]),]
      }
    }
    cat(x, "found", nrow(geneSNPCoupling[[x]]), "SNPs in gene\n")
    x <- x + 1
  }
  return(geneSNPCoupling)
}

BFMIsummary <- summarizeImprintedSNPsInGene(uniqueGenes, geneExonsCoupling, matBFMIsnps)
names(BFMIsummary) <- uniqueGenes
B6Nsummary  <- summarizeImprintedSNPsInGene(uniqueGenes, geneExonsCoupling, matB6Nsnps)
names(B6Nsummary) <- uniqueGenes

###
BFMIshort <- NULL
for(x in 1:length(BFMIsummary)){
  if(!is.null(BFMIsummary[[x]])) BFMIshort <- c(BFMIshort, BFMIsummary[x])
}

B6Nshort <- NULL
for(x in 1:length(B6Nsummary)){
  if(!is.null(B6Nsummary[[x]])) B6Nshort <- c(B6Nshort, B6Nsummary[x])
}

cat("Found",length(BFMIshort),"genes with usable",sum(unlist(lapply(BFMIshort,nrow))),"SNPs in matBFMI\n")         # Found 3681 genes with usable 20536 SNPs in matBFMI
cat("Found",length(B6Nshort), "genes with usable",sum(unlist(lapply(B6Nshort,nrow))), "SNPs in matB6N\n")          # Found 3503 genes with usable 19384 SNPs in matB6N

getASEgenes <- function(CROSSsummary, cutoff = 0.35) {
  v <- NULL
  for(x in 1:length(CROSSsummary)){ v <- c(v, mean(CROSSsummary[[x]][,"ImprintingScore"])); }
  hist(v)
  cat("Found", length(which(v > cutoff)),"ASE genes (cutoff =",cutoff,")\n")
  return(CROSSsummary[which(v > cutoff)])
}

B6Nase  <- getASEgenes(B6Nshort)                                                                                  # Found 59 ASE genes (cutoff = 0.35)
BFMIase <- getASEgenes(BFMIshort)                                                                                 # Found 93 ASE genes (cutoff = 0.35)

ASEmatrix <- NULL
for(x in 1:length(B6Nase)){
  ASEmatrix <- rbind(ASEmatrix, cbind(ensembl_gene_id = names(B6Nase[x]), cross = "matB6N", B6Nase[[x]]))
}
for(x in 1:length(BFMIase)){
  ASEmatrix <- rbind(ASEmatrix, cbind(ensembl_gene_id = names(BFMIase[x]), cross = "matBFMI", BFMIase[[x]]))
}

inRPKM <- match(ASEmatrix[,"ensembl_gene_id"],RPKM[,"ensembl_gene_id"])

refMean <- apply(ASEmatrix[,c("R1","R2","R3")],1,mean)
altMean <- apply(ASEmatrix[,c("A1","A2","A3")],1,mean)

ASEmatrix[,c("R1","R2","R3")] <- apply(ASEmatrix[,c("R1","R2","R3")], 2, function(x){round(x*100,2)})
ASEmatrix[,c("A1","A2","A3")] <- apply(ASEmatrix[,c("A1","A2","A3")], 2, function(x){round(x*100,2)})

ASEmatrix <- cbind(ensembl_gene_id = ASEmatrix[,"ensembl_gene_id"], cross =  ASEmatrix[,"cross"],  RPKM[inRPKM, c("mgi_symbol","mgi_description")], 
                   ASEmatrix[,-c(1:2)], refMean = round(refMean * 100, 2), altMean = round(altMean * 100, 2))

write.table(ASEmatrix, "ASE_10reads_noImputation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
matB6Nsnps  <- read.csv("maternalB6snps_10reads.txt", sep="\t", header=TRUE, colClasses="character")                                        # SNPs detected in the maternal B6N F1 cross
matBFMIsnps <- read.csv("maternalBFMIsnps_10reads.txt", sep="\t", header=TRUE, colClasses="character")                                      # SNPs detected in the maternal BFMI F1 cross
ASEmatrix <- read.table("ASE_10reads_noImputation.txt",sep="\t", header=TRUE, colClasses="character")

matB6Nfiles <- c("5070_CGATGT_L005_P_trimmed.aligned.sorted.realigned.recalibrated.bam", "5071_CCGTCC_L005_P_trimmed.aligned.sorted.realigned.recalibrated.bam", "5072_TAGCTT_L005_P_trimmed.aligned.sorted.realigned.recalibrated.bam")
matBFMIfiles <- c("5073_TTAGGC_L006_P_trimmed.aligned.sorted.realigned.recalibrated.bam", "5074_GATCAG_L006_P_trimmed.aligned.sorted.realigned.recalibrated.bam", "5075_ATGTCA_L006_P_trimmed.aligned.sorted.realigned.recalibrated.bam")

####

inMatBFMI <- unique(ASEmatrix[which(ASEmatrix[,"cross"] == "matBFMI"),"ID"])
inMatB6N <- unique(ASEmatrix[which(ASEmatrix[,"cross"] == "matB6N"),"ID"])

inBFMInotB6N <- inMatBFMI[which(!inMatBFMI %in% inMatB6N)]
inB6NnotBFMI <- inMatB6N[which(!inMatB6N %in% inMatBFMI)]

NoOtherSide <- c(as.character(inBFMInotB6N), as.character(inB6NnotBFMI))

recovered <- 0
imputed <- 0
notimputed <- 0

impM <- NULL
for(r in which(ASEmatrix[,"ID"] %in% NoOtherSide)){                                                   # Only the ones which are not duplicated need to be imputed
  chr <- ASEmatrix[r,"Chr"]
  pos <- ASEmatrix[r,"Loc"]
  ensembl_id <- as.character(ASEmatrix[r,"ensembl_gene_id"])
  if(ASEmatrix[r,"cross"] == "matBFMI"){
    inB6N <- which(as.character(matB6Nsnps[,"ID"]) == as.character(ASEmatrix[r,"ID"]))
    if(length(inB6N) > 0){                                                                            # No imputation, it was below the threshold
      impM <- rbind(impM, c(ensembl_id, "matB6N", ASEmatrix[r,"mgi_symbol"],ASEmatrix[r,"mgi_description"], matB6Nsnps[inB6N,], ASEmatrix[r,"Exon"], "", ""))
      recovered <- recovered + 1
    }else{                                                                                            # NEED TO IMPUTE, So look at the coverage
      coverage <- NULL
      for(infile in matB6Nfiles){
        cmd <- paste0("samtools mpileup -r ",chr,":",pos,"-",pos," Analysis/",infile)
        coverage <- c(coverage, as.numeric(strsplit(system(cmd, intern=TRUE)[1],"\t")[[1]][4]))
      }
      naCoverage <- which(is.na(coverage))
      if(length(naCoverage) > 0) coverage[naCoverage] <- 0
      info <- rep("", 12)
      if(sum(coverage > 10) == 3){ # Enough coverage
        info <- c(rep("BFMI", 3), rep("",9))
        if(ASEmatrix[r,"Detected"] == "BFMIsnp"){ info <- c(rep("B6N", 3), rep("",9)) }
        imputed <- imputed + 1
      }else{ notimputed <- notimputed + 1 }                                                            # Coverage too low
      impM <- rbind(impM, c(ensembl_id, "matB6N", ASEmatrix[r,"mgi_symbol"], ASEmatrix[r,"mgi_description"], 
      ASEmatrix[r,"ID"], ASEmatrix[r,"Chr"], ASEmatrix[r,"Loc"], ASEmatrix[r,"dbSNP"], info, ASEmatrix[r,"Detected"], ASEmatrix[r,"Exon"], "", ""))
      cat(coverage,"\n")
    }
  }
  if(ASEmatrix[r,"cross"] == "matB6N"){
    inBFMI <- which(as.character(matBFMIsnps[,"ID"]) == as.character(ASEmatrix[r,"ID"]))
    if(length(inBFMI) > 0){                                                                            # No imputation, it was below the threshold
      impM <- rbind(impM, c(ensembl_id, "matBFMI", ASEmatrix[r,"mgi_symbol"], ASEmatrix[r,"mgi_description"], matBFMIsnps[inBFMI,], ASEmatrix[r,"Exon"], "", ""))
      recovered <- recovered + 1
    }else{                                                                                             # NEED TO IMPUTE, So look at the coverage
      coverage <- NULL
      for(infile in matBFMIfiles){
        cmd <- paste0("samtools mpileup -r ",chr,":",pos,"-",pos," Analysis/",infile)
        coverage <- c(coverage, as.numeric(strsplit(system(cmd, intern=TRUE)[1],"\t")[[1]][4]))
      }
      naCoverage <- which(is.na(coverage))
      if(length(naCoverage) > 0) coverage[naCoverage] <- 0
      info <- rep("", 12)
      if(sum(coverage > 10) == 3){                                                                     # Enough coverage
        info <- c(rep("BFMI", 3), rep("",9))
        if(ASEmatrix[r,"Detected"] == "BFMIsnp"){ info <- c(rep("B6N", 3), rep("",9)) }
        imputed <- imputed + 1
      }else{ notimputed <- notimputed + 1 }                                                            # Coverage too low
      impM <- rbind(impM, c(ensembl_id, "matBFMI", ASEmatrix[r,"mgi_symbol"], ASEmatrix[r,"mgi_description"], 
      ASEmatrix[r,"ID"], ASEmatrix[r,"Chr"], ASEmatrix[r,"Loc"], ASEmatrix[r,"dbSNP"], info, ASEmatrix[r,"Detected"], ASEmatrix[r,"Exon"], "", ""))
      cat(coverage,"\n")
    }
  }
}
cat("We recovered", recovered, "SNPs, which were below the threshold\n")              # We recovered 24 SNPs, which were below the threshold
cat("We imputed", imputed, "SNPs since they had enough coverage\n")                   # We imputed 150 SNPs since they had enough coverage
cat("We failed", notimputed, "SNPs since they had not enough coverage\n")             # We failed 35 SNPs since they had not enough coverage
colnames(impM) <- colnames(ASEmatrix)

write.table(impM, "ASE_10reads_imputation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#
ASEmatrix <- read.table("ASE_10reads_noImputation.txt",sep="\t", header=TRUE)
IMPmatrix <- read.table("ASE_10reads_imputation.txt",sep="\t", header=TRUE)

refMean <- apply(IMPmatrix[,c("R1","R2","R3")],1,mean)
altMean <- apply(IMPmatrix[,c("A1","A2","A3")],1,mean)

IMPmatrix[,c("R1","R2","R3")] <- apply(IMPmatrix[,c("R1","R2","R3")],2,function(x){round(x*100,2)})
IMPmatrix[,c("A1","A2","A3")] <- apply(IMPmatrix[,c("A1","A2","A3")],2,function(x){round(x*100,2)})

IMPmatrix[,"refMean"] <- round(refMean * 100, 2)
IMPmatrix[,"altMean"] <- round(altMean * 100, 2)

ALLmatrix <- rbind(ASEmatrix,IMPmatrix)

write.table(ALLmatrix, "ASE_10reads_combined.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### Write the table out for the paper ###


ALLmatrixLongIDs <- paste0(ALLmatrix[,"ID"],"_",ALLmatrix[,"mgi_symbol"])
uniqueIDs <- unique(ALLmatrixLongIDs)
columns <- c("snpID", "dbSNP", "Chr", "Loc", "SNPorigin", "ensembl_gene_id", "mgi_symbol", "mgi_description", "exon", "ref_B6N", "alt_B6N", "score_B6N", "ASE_B6N", "ref_BFMI", "alt_BFMI", "score_BFMI", "ASE_BFMI", "Class")

paperMatrix <- matrix(NA, length(uniqueIDs), length(columns))
colnames(paperMatrix) <- columns

originalIDs <- unlist(lapply(lapply(strsplit(uniqueIDs,"_"),"[",c(1:2)),function(x){paste0(x[1],"_",x[2])}))

paperMatrix[,"snpID"] <- as.character(originalIDs)
staticLocs <- match(uniqueIDs, ALLmatrixLongIDs)
paperMatrix[,"dbSNP"]             <- as.character(ALLmatrix[staticLocs, "dbSNP"])
paperMatrix[,"Chr"]               <- as.character(ALLmatrix[staticLocs, "Chr"])
paperMatrix[,"Loc"]               <- ALLmatrix[staticLocs, "Loc"]
paperMatrix[,"SNPorigin"]         <- as.character(ALLmatrix[staticLocs, "Detected"])
paperMatrix[,"ensembl_gene_id"]   <- as.character(ALLmatrix[staticLocs, "ensembl_gene_id"])
paperMatrix[,"mgi_symbol"]        <- as.character(ALLmatrix[staticLocs, "mgi_symbol"])
paperMatrix[,"mgi_description"]   <- as.character(ALLmatrix[staticLocs, "mgi_description"])
paperMatrix[,"exon"]              <- as.character(ALLmatrix[staticLocs, "Exon"])


for(r in 1:nrow(paperMatrix)){
  ALLsubset <- ALLmatrix[which(ALLmatrix[,"ID"] == paperMatrix[r, "snpID"]),]
  matB6N <- which(ALLsubset[,"cross"]=="matB6N")
  paperMatrix[r,"ref_B6N"]   <- ALLsubset[matB6N, "refMean"][1]
  paperMatrix[r,"alt_B6N"]   <- ALLsubset[matB6N, "altMean"][1]
  paperMatrix[r,"score_B6N"] <- ALLsubset[matB6N, "ImprintingScore"][1]
  imprB6N <- "NotConsistent"
  if(ALLsubset[matB6N,"Origin1"] == ALLsubset[matB6N,"Origin2"] && ALLsubset[matB6N,"Origin1"] == ALLsubset[matB6N,"Origin3"]) imprB6N <- ALLsubset[matB6N,"Origin1"][1]
  paperMatrix[r,"ASE_B6N"] <- as.character(imprB6N)
  
  matBFMI <- which(ALLsubset[,"cross"]=="matBFMI")
  paperMatrix[r,"ref_BFMI"]   <- ALLsubset[matBFMI, "refMean"][1]
  paperMatrix[r,"alt_BFMI"]   <- ALLsubset[matBFMI, "altMean"][1]
  paperMatrix[r,"score_BFMI"] <- ALLsubset[matBFMI, "ImprintingScore"][1]
  imprBFMI <- "NotConsistent"
  if(ALLsubset[matBFMI,"Origin1"] == ALLsubset[matBFMI,"Origin2"] && ALLsubset[matBFMI,"Origin1"] == ALLsubset[matBFMI,"Origin3"]) imprBFMI <- ALLsubset[matBFMI,"Origin1"][1]
  paperMatrix[r,"ASE_BFMI"] <- as.character(imprBFMI)
  if(imprBFMI == "BFMI" && imprB6N == "BFMI") paperMatrix[r,"Class"] <- "BFMI"
  if(imprBFMI == "B6N" && imprB6N == "B6N") paperMatrix[r,"Class"] <- "B6N"
  if(imprBFMI == "B6N" && imprB6N == "BFMI") paperMatrix[r,"Class"] <- "Paternal"
  if(imprBFMI == "BFMI" && imprB6N == "B6N") paperMatrix[r,"Class"] <- "Maternal"
}

paperMatrix[is.na(paperMatrix)] <- ""

write.table(paperMatrix, "ASE_10reads_forPaper.txt", sep = "\t", quote = FALSE, row.names = FALSE)


################# OLD #######################
imputeReferenceASE <- function(BFMIsummary, B6Nsummary, RPKM, RPKMcutoff = 3, ASEcutoff = 0.35){
  for(x in 1:length(BFMIsummary)){
    if(!is.null(BFMIsummary[[x]]) && is.null(B6Nsummary[[x]])){
      if(mean(BFMIsummary[[x]][,"ImprintingScore"]) >= ASEcutoff && RPKM[which(RPKM[,"ensembl_gene_id"] == names(BFMIsummary[x])),"Mean.B6NxBFMI860.12.L"] >= RPKMcutoff){
        B6Nsummary[[x]] <- BFMIsummary[[x]]
        for(snp in 1:nrow(B6Nsummary[[x]])){
          for(origin in c("Origin1", "Origin2", "Origin3")){
            if(B6Nsummary[[x]][snp,"Detected"]=="BFMIsnp"){ B6Nsummary[[x]][snp,origin] <- "B6N"; }else{ B6Nsummary[[x]][snp,origin] <- "BFMI"; }
          }
          B6Nsummary[[x]][snp,"ImprintingScore"] <- 1
          for(column in c("R1", "A1", "R2", "A2", "R3", "A3")){ B6Nsummary[[x]][snp, column] <- "?"; }
        }
        cat(x,"imputed\n");
      }
    }
    if(is.null(BFMIsummary[[x]]) && !is.null(B6Nsummary[[x]])){
      if(mean(B6Nsummary[[x]][,"ImprintingScore"]) >= ASEcutoff && RPKM[which(RPKM[,"ensembl_gene_id"] == names(BFMIsummary[x])),"Mean.BFMI860.12xB6N.L"] >= RPKMcutoff){
        BFMIsummary[[x]] <- B6Nsummary[[x]]
        for(snp in 1:nrow(BFMIsummary[[x]])){
          for(origin in c("Origin1", "Origin2", "Origin3")){
            if(BFMIsummary[[x]][snp,"Detected"]=="BFMIsnp"){ BFMIsummary[[x]][snp,origin] <- "B6N"; }else{ BFMIsummary[[x]][snp,origin] <- "BFMI"; }
          }
          BFMIsummary[[x]][snp,"ImprintingScore"] <- 1
          for(column in c("R1", "A1", "R2", "A2", "R3", "A3")){ BFMIsummary[[x]][snp, column] <- "?"; }
        }
        cat(x,"imputed\n");
      }
    }
  }
  return(list(BFMIsummary, B6Nsummary))
}

imputedData <- imputeReferenceASE(BFMIsummary, B6Nsummary, RPKM)                                                 # Impute the expressed genes (without a SNP)


BFMIase <- getShortList(imputedData[[1]])
for(x in 1:length(BFMIase)){ write.table(BFMIase[[x]], paste0(names(BFMIase[x]),".BFMI.txt"),sep="\t", quote=FALSE); }

B6Nase <- getShortList(imputedData[[2]])
for(x in 1:length(B6Nase)){ write.table(B6Nase[[x]], paste0(names(B6Nase[x]),".B6N.txt"),sep="\t", quote=FALSE); }


summarize <- function(CROSSsummary){
  mmatrix <- matrix(NA, length(CROSSsummary), 4)                                                                # Create an output matrix with 3 columns
  for(x in 1:length(CROSSsummary)){ 
    mmatrix[x,1] <- names(CROSSsummary[x])                                                                      # Ensembl geneID
    mmatrix[x,2] <- mean(CROSSsummary[[x]][,"ImprintingScore"]);                                                # Mean imprinting score across the gene
    origin <- names(which.max(table(unlist(CROSSsummary[[x]][,c("Origin1","Origin2","Origin3")]))))
    if(!is.null(origin)){ mmatrix[x,3] <- origin; }                                                             # Take the one which occurs most as the origin
    detected <- names(which.max(table(unlist(CROSSsummary[[x]][,c("Detected")]))))
    if(!is.null(origin)){ mmatrix[x,4] <- detected; }                                                           # Take the one which occurs most as the origin
  }
  colnames(mmatrix) <- c("ensembl_gene_id", "ImprintingScore", "Origin", "Detected")
  return(mmatrix)
}

BFMIaseSummary <- summarize(BFMIase)
B6NaseSummary <- summarize(B6Nase)

# TODO: Filter the shortList for possible errors due to the other side not expressing
filterASE <- function(CROSSsummary, otherCross, RPKM, RPKMcutoff = 3, noFilter=FALSE){
  ordering <- match(CROSSsummary[,"ensembl_gene_id"], RPKM[,"ensembl_gene_id"])
  combined <- cbind(RPKM[ordering,], CROSSsummary)
  if(noFilter)return(combined)
  cat("Starting with", nrow(combined),"\n")
  ## BFMI based genes, we should remove because of low expression in B6N
  SnpFromBFMInoB6Nexpression <-  which(combined[,"Origin"] == "BFMI" & as.numeric(combined[,"Mean.B6N"]) < RPKMcutoff)
  genesWeCouldRemove <- combined[SnpFromBFMInoB6Nexpression, "ensembl_gene_id"]
  genesToRemove <- NULL
  for(gene in genesWeCouldRemove){
    imp1 <- CROSSsummary[which(CROSSsummary[,"ensembl_gene_id"] == gene), "Origin"]                                               # Only delete genes which have the same imprinted SNP
    imp2 <- otherCross[which(otherCross[,"ensembl_gene_id"] == gene), "Origin"]
    if(length(imp2) > 0 && imp1 == imp2) genesToRemove <- c(genesToRemove, gene)
  }
  combined <- combined[ - which(combined[,"ensembl_gene_id"] %in% genesToRemove),]                                                # Filter the erroneously called BFMI ASE genes
  
  ## B6N based genes, we should remove because of low expression in BFMI
  SnpFromBFMInoB6Nexpression <-  which(combined[,"Origin"] == "B6N" & as.numeric(combined[,"Mean.BFMI860"]) < RPKMcutoff)
  genesWeCouldRemove <- combined[SnpFromBFMInoB6Nexpression, "ensembl_gene_id"]
  genesToRemove <- NULL
  for(gene in genesWeCouldRemove){
    imp1 <- CROSSsummary[which(CROSSsummary[,"ensembl_gene_id"] == gene), "Origin"]                                               # Only delete genes which have the same imprinted SNP
    imp2 <- otherCross[which(otherCross[,"ensembl_gene_id"] == gene), "Origin"]
    if(length(imp2) > 0 && imp1 == imp2) genesToRemove <- c(genesToRemove, gene)
  }
  combined <- combined[ - which(combined[,"ensembl_gene_id"] %in% genesToRemove),]                                                # Filter the erroneously called B6N ASE genes
  
  cat("Filtered out", nrow(CROSSsummary) - nrow(combined), "genes due to low expression in one of the parents\n")
  predictedGenes <- grep("predicted", combined[,"mgi_description"])
  combined <- combined[ - predictedGenes, ]
  cat("Filtered out", length(predictedGenes), "predicted genes, left with",  nrow(combined),"\n")
  return(combined)
}

write.table(filterASE(BFMIaseSummary, B6NaseSummary, RPKM), "ASE_matBFMIsnps_5reads.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(filterASE(B6NaseSummary, BFMIaseSummary, RPKM), "ASE_matB6Nsnps_5reads.txt", sep="\t", quote=FALSE, row.names=FALSE)

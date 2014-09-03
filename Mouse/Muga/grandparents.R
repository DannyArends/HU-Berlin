# Phasing the MUGA SNP array towards grandparents
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sep, 2014
# first written Sep, 2014

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                  # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

parentinfo <- phenotypedata[,c("Vater", "Mutter")]

setwd("E:/Mouse/DNA/MegaMuga/")
genotypes <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE, colClasses="character")          # Genotypes measured on the MUGA array

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })                            # Missing amount of genotype data per individual
genotypes <- genotypes[,-which(missingPerInd==100)]                                                             # Remove individuals which have NO genotypes

parentinfo <- parentinfo[which(rownames(parentinfo) %in% colnames(genotypes)),]

individual <- "6661051"

hasParents <- function(parentinfo, individual = "6661051"){
  if(any(parentinfo[individual, ] == 0)){ return(FALSE) }                                                       # We do not have parents in the set
  return(TRUE)
}

isGrandChild <- function(parentinfo, individual = "6661051"){
  if(hasParents(parentinfo, individual)){
    for(parent in parentinfo[individual, ]){ if(!hasParents(parentinfo, as.character(parent))){ return(FALSE); } }
    return(TRUE)
  }
  return(FALSE)
}

getGP <- function(parentinfo, individual = "6661051", side="Vater"){                                            # Get the grandparents *and parent as rowname*
  if(!isGrandChild(parentinfo, individual)) stop(paste0("individual", individual, "has no grandparents"))
  return(parentinfo[as.character(parentinfo[individual, side]),])
}

phase <- function(genotypes, parentinfo, individual, verbose = TRUE){
  if(!isGrandChild(parentinfo, individual)) stop(paste0("individual", individual, "has no grandparents"))
  paternal <- getGP(parentinfo, individual, "Vater")
  maternal <- getGP(parentinfo, individual, "Mutter")
  mfamily <- unlist(c(individual, rownames(paternal), paternal, rownames(maternal), maternal))
  phasedgenotype <- rep("", nrow(genotypes))
  for(m in 1:nrow(genotypes)){
    if(!any(is.na(genotypes[m,mfamily]))){
      if(genotypes[m,individual] == "H"){
        #cat(m, "Individual is",genotypes[m,individual],", Dad =", genotypes[m, rownames(paternal)], ", Mom =", genotypes[m, rownames(maternal)],"\n")
        if(genotypes[m, rownames(paternal)] == "H" && genotypes[m, rownames(maternal)] != "H"){
          GF <- as.character(paternal[,"Vater"])
          GM <- as.character(paternal[,"Mutter"])
          if(genotypes[m, GF] == "A" && genotypes[m, GM] != "A"){ if(verbose){ cat(m, "Allele from: GrandMother of Fathers side\n"); }; phasedgenotype[m] <- "GmF"; }
          if(genotypes[m, GF] == "B" && genotypes[m, GM] != "B"){ if(verbose){ cat(m, "Allele from: GrandMother of Fathers side\n"); }; phasedgenotype[m] <- "GmF"; }
          if(genotypes[m, GM] == "A" && genotypes[m, GF] != "A"){ if(verbose){ cat(m, "Allele from: GrandFather of Fathers side\n"); }; phasedgenotype[m] <- "GfF"; }
          if(genotypes[m, GM] == "B" && genotypes[m, GF] != "B"){ if(verbose){ cat(m, "Allele from: Grandfather of Fathers side\n"); }; phasedgenotype[m] <- "GfF"; }
        }
        if(genotypes[m, rownames(maternal)] == "H" && genotypes[m, rownames(paternal)] != "H"){
          GF <- as.character(maternal[,"Vater"])
          GM <- as.character(maternal[,"Mutter"])
          if(genotypes[m, GF] == "A" && genotypes[m, GM] != "A"){ if(verbose){ cat(m, "Allele from: GrandMother of Mothers side\n"); }; phasedgenotype[m] <- "GmM"; }
          if(genotypes[m, GF] == "B" && genotypes[m, GM] != "B"){ if(verbose){ cat(m, "Allele from: GrandMother of Mothers side\n"); }; phasedgenotype[m] <- "GmM"; }
          if(genotypes[m, GM] == "A" && genotypes[m, GF] != "A"){ if(verbose){ cat(m, "Allele from: GrandFather of Mothers side\n"); }; phasedgenotype[m] <- "GfM"; }
          if(genotypes[m, GM] == "B" && genotypes[m, GF] != "B"){ if(verbose){ cat(m, "Allele from: Grandfather of Mothers side\n"); }; phasedgenotype[m] <- "GfM"; }
        }
      }
    }
  }
  return(phasedgenotype)
}

phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),]       # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]

phasedgenotypes <- NULL
for(individual in F2){
  phasedgenotypes <- cbind(phasedgenotypes, phase(genotypes, parentinfo, individual, FALSE))
  cat("Done individual", individual, "\n")
}

colnames(phasedgenotypes) <- F2
rownames(phasedgenotypes) <- rownames(genotypes)
write.table(phasedgenotypes, "Analysis/genotypesPhasedGP.txt", sep="\t")

phasedSubSet <- phasedgenotypes[which(apply(phasedgenotypes,1,function(x){sum(x!="")}) > 50),]
write.table(phasedSubSet, "Analysis/genotypesPhasedGPSubSet.txt", sep="\t")

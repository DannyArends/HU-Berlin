# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

library("psych")

source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes
genotypesPh <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("", "AA", "CC", "TT", "GG", "??"))    # Phased by Beagle (only heterozygous)
genotypesFP <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("", "??"))                            # Phased by Beagle (ALL)
genotypesGP <- read.table("Analysis/genotypesPhasedGP.txt",     sep="\t", check.names=FALSE, colClasses="character", na.strings=c(""))                                  # Phased towards the grandparents

phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals

corG <- read.table("Analysis/genotypecorrelation.txt", sep = "\t", check.names=FALSE)                         # Genotype correlation matrix
pcas <- principal(corG, nfactors = 1)                                                                         # Calculate the first principal component

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
birthmonth <- unlist(lapply(strsplit(as.character(phenotypes[,"W.dat"]),".", fixed=TRUE),"[",2))
phenotypes <- cbind(phenotypes, Birthmonth = birthmonth)                                                      # Add the birth month column to the matrix

phenotypes <- cbind(phenotypes, mri42d_fatDlean = phenotypes[,"mri42d_fat"] / phenotypes[,"mri42d_lean"])     # Fat / Lean day 42
phenotypes <- cbind(phenotypes, mri56d_fatDlean = phenotypes[,"mri56d_fat"] / phenotypes[,"mri56d_lean"])     # Fat / Lean day 56
phenotypes <- cbind(phenotypes, mri70d_fatDlean = phenotypes[,"mri70d_fat"] / phenotypes[,"mri70d_lean"])     # Fat / Lean day 70
phenotypes <- cbind(phenotypes, PC1 = NA)                                                                     # Empty PC1
phenotypes[names(pcas$loadings[,1]), "PC1"] <- pcas$loadings[,1]

phenotypes <- cbind(phenotypes, mri70d_fatDlean = phenotypes[,"mri70d_fat"] / phenotypes[,"mri70d_lean"])     # Fat / Lean day 70

onegenotype <- which(lapply(apply(genotypes[,F2], 1, table), length) == 1)                                    # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- unique(genotypes[-onegenotype, F2])                                                            # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypes), "markers\n")                                                                # == Left with 11581 markers

onegenotype <- which(lapply(apply(genotypesPh[,F2], 1, table), length) == 1)                                  # Markers with only one genotype cannot be used in QTL mapping
genotypesPh <- unique(genotypesPh[-onegenotype,F2])                                                           # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesPh), "markers\n")                                                              # == Left with 19465 markers

onegenotype <- which(lapply(apply(genotypesFP[,F2], 1, table), length) == 1)                                  # Markers with only one genotype cannot be used in QTL mapping
genotypesFP <- unique(genotypesFP[-onegenotype,F2])                                                           # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesFP), "markers\n")                                                              # == Left with 19908 markers

onegenotype <- which(lapply(apply(genotypesGP[,F2], 1, table), length) == 1)                                  # Markers with only one genotype cannot be used in QTL mapping
genotypesGP <- unique(genotypesGP[-onegenotype,F2])                                                           # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesGP), "markers\n")                                                              # == Left with 7551 markers

mriGWAS <- function(genotypes, phenotypes, pheno.col = "42d", to = nrow(genotypes), cof = ""){                # to = nrow(genotypes)
  nterms  <- 5;
  if(cof != "") nterms <- 6
  pvalues <- matrix(NA, min(to, nrow(genotypes)), nterms)
  rownames(pvalues) <- rownames(genotypes)[1:to]
  if(cof == ""){
    colnames(pvalues) <- c("subfamily", "l_size", "l_number", "season", "marker")
  }else{
    colnames(pvalues) <- c("subfamily", "l_size", "l_number", "season", cof, "marker")
  }
  for(x in 1:to){
    ind            <- colnames(genotypes[x,!is.na(genotypes[x,])])                                                            # Which individuals have genotype data

    subfamily      <- as.numeric(phenotypes[ind, "PC1"])                                                                     # Fixed effect: Subfamily structure (factor)
    littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter  (linear effect)
    litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter    (factor)
    season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born    (factor)
    genotype       <- as.factor(t(genotypes[x,!is.na(genotypes[x,])]))                                                        # The genotype under investigation  (factor)

    phenotype      <- phenotypes[ind, pheno.col]                      #/ phenotypes[ind, paste0("mri",pheno.col,"_lean")]     # Response: Fat / Lean
    if(cof == ""){
      tryCatch(
        res <- anova(lm(phenotype ~ subfamily + littersize + litternumber + season + genotype))[[5]], 
        error = function(e){ res <<- rep(NA, 5) })
    }else{
      covar        <- as.factor(t(genotypes[cof, ind]))
      tryCatch(res <- anova(lm(phenotype ~ subfamily + littersize + litternumber + season + covar + genotype))[[5]], error = function(e){ res <<- rep(NA, 6) })
    }
    if(x %% 1000 == 0) cat("Marker:", x, round(-log10(res[-length(res)]),1),"\n")
    pvalues[x, 1:(length(res)-1)] <- res[-length(res)]
  }
  return(round(-log10(pvalues), 3))
}

phenonames <- c("mri42d_fatDlean", "mri56d_fatDlean", "mri70d_fatDlean",                                        # Fat / Lean 
                "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71",                                  # Weights
                "GF1", "GF2", "total.GF",                                                                       # Gonadal fat
                "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",                                 # Other phenotypes
                "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4",                                            # MRI day 42
                "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4",                                            # MRI day 56
                "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4")                                            # MRI day 70

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
for(phe in phenonames){
  qtls <- mriGWAS(genotypes,   phenotypes, phe);                                                                # Map QTLs, normal GWAS (A, H, B)
  write.table(qtls, paste0("Analysis/PCA_qtls_", phe, ".txt"),   sep="\t")                                          # Write results
  cat("QTLs done for", phe, "\n")

#  qtlsC   <- mriGWAS(genotypes,   phenotypes, phe, cof = "UNC5048297") ;  
#  write.table(qtlsC, paste0("Analysis/PCA_qtls_", phe, "_cof_UNC5048297.txt"),   sep="\t")                          # Write results
#  cat("QTL + COF done for", phe, "\n")
  
#  phased_qtls <- mriGWAS(genotypesPh, phenotypes, phe) ;                                                        # Map QTLs, GWAS on the H0 versus H1
#  write.table(phased_qtls, paste0("Analysis/PCA_qtls_phased_", phe, ".txt"),   sep="\t")                            # Write results
#  cat("QTL PHASED done for", phe, "\n")

#  phasedfull_qtls <- mriGWAS(genotypesFP, phenotypes, phe) ;                                                    # Map QTLs, GWAS on full phase: A, H0, H1, B
#  write.table(phasedfull_qtls, paste0("Analysis/PCA_qtls_phasedfull_", phe, ".txt"),   sep="\t")                    # Write results
#  cat("QTL PHASED FULL done for", phe, "\n")

#  phasedfull_qtlsC <- mriGWAS(genotypes,   phenotypes, phe, cof = "UNC5048297") ; 
#  write.table(phasedfull_qtlsC, paste0("Analysis/PCA_qtls_phasedfull_", phe, "_cof_UNC5048297.txt"),   sep="\t")    # Write results
#  cat("QTL PHASED FULL + COF done for", phe, "\n")
  
#  grandparents_qtls <- mriGWAS(genotypesGP, phenotypes, phe) ;                                                  # Map QTLs, GWAS towards the grandparents
#  write.table(phased_qtls, paste0("Analysis/PCA_qtls_grandparents_", phe, ".txt"),   sep="\t")                      # Write results
#  cat("QTL PHASED GP done for", phe, "\n")
}

plotZoom <- function(qtls, smap, chr){
  onChr <- rownames(smap[which(smap[,"Chr"] == chr),])
  plot(as.numeric(smap[onChr,"Mb_NCBI38"]) / 1000000, qtls[onChr,"BH"], t = 'p', xlab = paste0("Chromosome", chr), main = paste0("QTL profile ", phe, " (zoom)"), ylab="LOD", pch=19, cex=0.2,col="blue")
  ma <- function(x, n = 5){ filter(x, rep(1/n,n), sides=2) }
  points(as.numeric(smap[onChr,"Mb_NCBI38"]) / 1000000, ma(qtls[onChr,"BH"]), t = 'l', main = paste0("QTL profile, (zoom Chr ",chr,")"), ylab="LOD",lwd=3)
  points(as.numeric(smap[onChr,"Mb_NCBI38"]) / 1000000, qtls[onChr,"BH"], t = 'o', main = paste0("QTL profile, (zoom Chr ",chr,")"), ylab="LOD", pch=19, cex=0.2, col="blue",lty=2)
  abline(h = 3, col = "orange", lty=2) ; abline(h = 5, col = "green", lty=2)
}


setwd("E:/Mouse/DNA/MegaMuga/")
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
chrcolors <- rep(c("black","orange"),length(unique(map[,"Chr"])))
names(chrcolors) <- unique(map[,"Chr"])

setwd("E:/Mouse/ClassicalPhenotypes/AIL")

for(phe in phenonames){
  qtls <- read.table(paste0("Analysis/PCA_qtls_", phe, ".txt"),sep="\t")
  qtls <- cbind(qtls, BH  = round(-log10(p.adjust(10^(-qtls[,"marker"]), "BH")), d = 2))

  significant <- rownames(qtls[which(qtls[,"BH"] > 3),])
  if(length(significant) > 0){
    cat("---", phe,"\n")
    smap <- map[rownames(qtls),]
    cat("QTLs on chromosome:", unique(smap[significant,"Chr"]),"\n")
    plot(qtls[,"BH"], t='h', col=chrcolors[smap[rownames(qtls),"Chr"]], main = paste("QTL profile", phe), ylab="LOD")
    
    for(chr in unique(smap[significant,"Chr"])){
      cat("Chromosome", chr, "\n")
      onChr <- rownames(smap[which(smap[,"Chr"] == chr),])
      top <- which.max(qtls[onChr,"BH"])
      for(x in max(1,(top-10)):(top+10)){
        mName <- rownames(qtls[onChr[x],])
        cat(x,mName, smap[mName,"Mb_NCBI38"], qtls[mName,"marker"], qtls[mName,"BH"],"\n")
        #cat(x, mName,"\n")
      }
      plotZoom(qtls, smap, chr)
    }
  }else{
    cat("---", phe,"no significant QTL\n")
  }
}

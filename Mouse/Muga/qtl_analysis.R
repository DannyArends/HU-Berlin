# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2015
# first written Aug, 2014

source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_paper/dateToSeason.R")
setwd("D:/Edrive//Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes
genotypesPh <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("", "AA", "CC", "TT", "GG", "??"))    # Phased by Beagle (only heterozygous)
genotypesFP <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("", "??"))                            # Phased by Beagle (ALL)
genotypesGP <- read.table("Analysis/genotypesPhasedGP.txt",     sep="\t", check.names=FALSE, colClasses="character", na.strings=c(""))                                  # Phased towards the grandparents

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("33310233", "6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                         # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]; length(F2)                                     # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]; length(F1)                                     # The F1 individuals
P <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]; length(P)                                       # The P individuals

## Lecture
#subsetMkr <- sample(rownames(subsetG),3000)

#subsetP <- phenotypes[F2, c(paste0("d",seq(21,70,7)),"WG2")]
#subsetG <- genotypes[subsetMkr, F2]
#subsetM <- map[subsetMkr,]
#setwd("D:/Projects/Lectures 2014-2015/R course/R course/Assignments")

#write.table(subsetP, "lecture5_phenotypes.txt", sep="\t", quote=FALSE)
#write.table(subsetG, "lecture5_genotypes.txt", sep="\t", quote=FALSE)
#write.table(subsetM, "lecture5_map.txt", sep="\t", quote=FALSE)

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
birthmonth <- unlist(lapply(strsplit(as.character(phenotypes[,"W.dat"]),".", fixed=TRUE),"[",2))
phenotypes <- cbind(phenotypes, Birthmonth = birthmonth)                                                      # Add the birth month column to the matrix

phenotypes <- cbind(phenotypes, mri42d_fatDlean = phenotypes[,"mri42d_fat"] / phenotypes[,"mri42d_lean"])     # Fat / Lean day 42
phenotypes <- cbind(phenotypes, mri56d_fatDlean = phenotypes[,"mri56d_fat"] / phenotypes[,"mri56d_lean"])     # Fat / Lean day 56
phenotypes <- cbind(phenotypes, mri70d_fatDlean = phenotypes[,"mri70d_fat"] / phenotypes[,"mri70d_lean"])     # Fat / Lean day 70

onegenotype <- which(lapply(apply(genotypes[,F2], 1, table), length) < 2)                                       # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- (genotypes[-onegenotype, F2])                                                                    # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypes), "markers\n")                                                                  # == Left with 21260 markers    // Left with 11581 markers

onegenotype <- which(lapply(apply(genotypesPh[,F2], 1, table), length) < 2)                                     # Markers with only one genotype cannot be used in QTL mapping
genotypesPh <- (genotypesPh[-onegenotype,F2])                                                                   # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesPh), "markers\n")                                                                # == Left with 20199 markers    // Left with 19465 markers

onegenotype <- which(lapply(apply(genotypesFP[,F2], 1, table), length) < 2)                                     # Markers with only one genotype cannot be used in QTL mapping
genotypesFP <- (genotypesFP[-onegenotype,F2])                                                                   # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesFP), "markers\n")                                                                # == Left with 20930 markers    // Left with 19908 markers

onegenotype <- which(lapply(apply(genotypesGP[,F2], 1, table), length) < 2)                                     # Markers with only one genotype cannot be used in QTL mapping
genotypesGP <- (genotypesGP[-onegenotype,F2])                                                                   # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesGP), "markers\n")                                                                # == Left with 19903 markers    // Left with 7551 markers

phenotypes <- phenotypes[F2,]

#Growth curves
growth <- phenotypes[F2, paste0("d",seq(21,70,7))]
plot(c(1,length(seq(21,70,7))), c(0,max(growth)), t = 'n', main="Growth curves AIL individuals", sub="Generation 28", ylab="Bodyweight (g)", xlab="Day",xaxt='n')
colorz <- c("orange", "black", "gray")[as.numeric(as.factor(unlist(genotypes["UNC5048297", F2])))]
x <- 1
apply(growth,1,function(x){points(x, t = 'l', col=colorz[x]); x <<- x + 1 } )


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

    subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure (factor)
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

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
for(phe in phenonames){
#  qtls <- mriGWAS(genotypes,   phenotypes, phe);                                                                # Map QTLs, normal GWAS (A, H, B)
#  write.table(qtls, paste0("Analysis/qtls_", phe, ".txt"),   sep="\t")                                          # Write results
#  cat("QTLs done for", phe, "\n")

#  qtlsC   <- mriGWAS(genotypes,   phenotypes, phe, cof = "UNC5048297") ;  
#  write.table(qtlsC, paste0("Analysis/qtls_", phe, "_cof_UNC5048297.txt"),   sep="\t")                          # Write results
#  cat("QTL + COF done for", phe, "\n")
  
#  phased_qtls <- mriGWAS(genotypesPh, phenotypes, phe) ;                                                        # Map QTLs, GWAS on the H0 versus H1
#  write.table(phased_qtls, paste0("Analysis/qtls_phased_", phe, ".txt"),   sep="\t")                            # Write results
#  cat("QTL PHASED done for", phe, "\n")

#  phasedfull_qtls <- mriGWAS(genotypesFP, phenotypes, phe) ;                                                    # Map QTLs, GWAS on full phase: A, H0, H1, B
#  write.table(phasedfull_qtls, paste0("Analysis/qtls_phasedfull_", phe, ".txt"),   sep="\t")                    # Write results
#  cat("QTL PHASED FULL done for", phe, "\n")

#  phasedfull_qtlsC <- mriGWAS(genotypes,   phenotypes, phe, cof = "UNC5048297") ; 
#  write.table(phasedfull_qtlsC, paste0("Analysis/qtls_phasedfull_", phe, "_cof_UNC5048297.txt"),   sep="\t")    # Write results
#  cat("QTL PHASED FULL + COF done for", phe, "\n")
  
#  grandparents_qtls <- mriGWAS(genotypesGP, phenotypes, phe) ;                                                  # Map QTLs, GWAS towards the grandparents
#  write.table(grandparents_qtls, paste0("Analysis/qtls_grandparents_", phe, ".txt"),   sep="\t")                # Write results
#  cat("QTL PHASED GP done for", phe, "\n")
}

#for(phe in phenonames){
#  phasedFULLC <- read.table(paste0("Analysis/qtls_phasedfull_", phe, "_cof_UNC5048297.txt"),   sep="\t")
#  plot(phasedFULLC[,"marker"], main=phe)
#  scan()
#}

plotZoom <- function(qtls, smap, chr){
  onChr <- rownames(smap[which(smap[,"Chr"] == chr),])[1:500]
  plot(as.numeric(smap[onChr,"Mb_NCBI38"]) / 1000000, qtls[onChr,"BH"], t = 'p', xlab = paste0("Chromosome", chr), main = paste0("QTL profile ", phe, " (zoom)"), ylab="LOD", pch=19, cex=0.2,col="blue")
  ma <- function(x, n = 5){ filter(x, rep(1/n,n), sides=2) }
  points(as.numeric(smap[onChr,"Mb_NCBI38"]) / 1000000, ma(qtls[onChr,"BH"]), t = 'l', main = paste0("QTL profile, (zoom Chr ",chr,")"), ylab="LOD",lwd=3)
  points(as.numeric(smap[onChr,"Mb_NCBI38"]) / 1000000, qtls[onChr,"BH"], t = 'o', main = paste0("QTL profile, (zoom Chr ",chr,")"), ylab="LOD", pch=19, cex=0.2, col="blue",lty=2)
  abline(h = 3, col = "orange", lty=2) ; abline(h = 5, col = "green", lty=2)
}

chrcolors <- rep(c("black","orange"),length(unique(map[,"Chr"])))
names(chrcolors) <- unique(map[,"Chr"])

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
for(phe in phenonames){
    qtls <- read.table(paste0("Analysis/qtls_", phe, ".txt"),sep="\t")
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
      for(x in max(1,(top-20)):(top+20)){
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

cbind(map[rownames(qtls[2856:2876,]),],qtls[2856:2876,])        # Chr 3 - 14 mb
cbind(map[rownames(qtls[13951:14131,]),],qtls[13951:14131,])    # Chr 12 - 9.7 mb

# ix <- which(rownames(qtls) == "UNC4784019")
# qtls[(ix-50):(ix+50),]



library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
bm.above <- getBM(attributes = c("ensembl_gene_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), 
filters = c("chromosomal_region", "biotype"), values = list("3:14841429:17006652","protein_coding"), mart = bio.mart)

bm.above <- getBM(attributes = c("ensembl_gene_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), 
filters = c("chromosomal_region", "biotype"), values = list("3:36599358:36762933","protein_coding"), mart = bio.mart)

bm.above <- getBM(attributes = c("ensembl_gene_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), 
filters = c("chromosomal_region", "biotype"), values = list("12:89085109:98002736","protein_coding"), mart = bio.mart)
write.table(bm.above, file="analysis/GenesinRegion_chr12.txt", sep="\t", row.names=FALSE)

bm.above <- getBM(attributes = c("ensembl_gene_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), 
filters = c("chromosomal_region", "biotype"), values = list("16:85457815:85919430","protein_coding"), mart = bio.mart)
write.table(bm.above, file="analysis/GenesinRegion_chr16.txt", sep="\t", row.names=FALSE)


qtl56   <- mriGWAS(genotypes,   phenotypes, "56d") ; 
qtl70   <- mriGWAS(genotypes,   phenotypes, "70d")  # Normal GWAS (A, H, B)
qtlPH42 <- mriGWAS(genotypesPh, phenotypes, "42d") ; 
qtlPH56 <- mriGWAS(genotypesPh, phenotypes, "56d") ; 
qtlPH70 <- mriGWAS(genotypesPh, phenotypes, "70d")  # GWAS on the H0 versus H1

qtl42C   <- mriGWAS(genotypes,   phenotypes, "42d", cof = "UNC5048297") ;
qtl56C   <- mriGWAS(genotypes,   phenotypes, "56d", cof = "UNC5048297") ;
qtl70C   <- mriGWAS(genotypes,   phenotypes, "70d", cof = "UNC5048297") ;


write.table(qtl42,   "Analysis/qtls_fatDlean42.txt",   sep="\t")

write.table(qtl42,   "Analysis/qtls_lean42.txt",   sep="\t")
write.table(qtl56,   "Analysis/qtls_fatDlean56.txt",   sep="\t")
write.table(qtl56,   "Analysis/qtls_fat56.txt",   sep="\t")
write.table(qtl56,   "Analysis/qtls_lean56.txt",   sep="\t")
write.table(qtl70,   "Analysis/qtls_fatDlean70.txt",   sep="\t")

write.table(qtl42C,   "Analysis/qtls_fatDlean42_UNC5048297.txt",   sep="\t")
write.table(qtl56C,   "Analysis/qtls_fatDlean56_UNC5048297.txt",   sep="\t")
write.table(qtl70C,   "Analysis/qtls_fatDlean70_UNC5048297.txt",   sep="\t")

write.table(qtlPH42, "Analysis/qtls_fatDlean42_PH.txt", sep="\t")
write.table(qtlPH56, "Analysis/qtls_fatDlean56_PH.txt", sep="\t")
write.table(qtlPH70, "Analysis/qtls_fatDlean70_PH.txt", sep="\t")

noChr3QTL <- which(as.character(genotypes["UNC5048297",F2]) != "A")

qtl42NoChr3IND   <- mriGWAS(genotypes[,F2[noChr3QTL]],   phenotypes[F2[noChr3QTL],], "42d") ;
write.table(qtl42NoChr3IND, "Analysis/qtls_fatDlean42_NoChr3IND.txt", sep="\t")

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
qtl42 <- read.table("Analysis/qtls_fatDlean42.txt",   sep="\t")
qtl56 <- read.table("Analysis/qtls_fatDlean56.txt",   sep="\t")
qtl70 <- read.table("Analysis/qtls_fatDlean70.txt",   sep="\t")

chrcolors <- rep(c("black","orange"),length(unique(map[,"Chr"])))
names(chrcolors) <- unique(map[,"Chr"])

png("Plots/QTL_FatdLean_Day42.png", width=1024, height=768)
  plot(qtl42[,"marker"], t='h', col=chrcolors[map[rownames(qtl42),"Chr"]], main = "QTL profile Fat/Lean day 42", ylab="LOD")
  abline(h = -log10(0.1/nrow(qtl42)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtl42)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtl42)), col="green", lty=2)
dev.off()

png("Plots/QTL_FatdLean_Day56.png", width=1024, height=768)
  plot(qtl56[,"marker"], t='h', col=chrcolors[map[rownames(qtl56),"Chr"]], main = "QTL profile Fat/Lean day 56", ylab="LOD")
  abline(h = -log10(0.1/nrow(qtl56)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtl56)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtl56)), col="green", lty=2)
dev.off()

png("Plots/QTL_FatdLean_Day70.png", width=1024, height=768)  
  plot(qtl70[,"marker"], t='h', col=chrcolors[map[rownames(qtl70),"Chr"]], main = "QTL profile Fat/Lean day 70", ylab="LOD")
  abline(h = -log10(0.1/nrow(qtl70)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtl70)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtl70)), col="green", lty=2)
dev.off()

qtl42C <- read.table("Analysis/qtls_fatDlean42_UNC5048297.txt",   sep="\t")
qtl56C <- read.table("Analysis/qtls_fatDlean56_UNC5048297.txt",   sep="\t")
qtl70C <- read.table("Analysis/qtls_fatDlean70_UNC5048297.txt",   sep="\t")

png("Plots/QTL_FatdLean_Day42_COF.png", width=1024, height=768)
  plot(qtl42C[,"marker"], t='h', col=chrcolors[map[rownames(qtl42C),"Chr"]], main = "QTL profile Fat/Lean day 42 (Cofactor)", ylab="LOD")
  abline(h = -log10(0.1/nrow(qtl42C)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtl42C)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtl42C)), col="green", lty=2)
dev.off()

png("Plots/QTL_FatdLean_Day56_COF.png", width=1024, height=768)
  plot(qtl56C[,"marker"], t='h', col=chrcolors[map[rownames(qtl56C),"Chr"]], main = "QTL profile Fat/Lean day 56 (Cofactor)", ylab="LOD")
  abline(h = -log10(0.1/nrow(qtl56C)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtl56C)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtl56C)), col="green", lty=2)
dev.off()

png("Plots/QTL_FatdLean_Day70_COF.png", width=1024, height=768)  
  plot(qtl70C[,"marker"], t='h', col=chrcolors[map[rownames(qtl70C),"Chr"]], main = "QTL profile Fat/Lean day 70 (Cofactor)", ylab="LOD")
  abline(h = -log10(0.1/nrow(qtl70C)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtl70C)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtl70C)), col="green", lty=2)
dev.off()


######


image(x = 1:nrow(qtls), y=(1:3)-0.5, as.matrix(qtls), oldstyle=TRUE, breaks=c(0,3,6,12,30,100), col=c("white",gray.colors(4)[4:1]))
grid(3); box()

map[which(qtls[,"qtl42"] > -log10(0.01/nrow(qtls))),]
map[which(qtls[,"qtl56"] > -log10(0.01/nrow(qtls))),]
map[which(qtls[,"qtl70"] > -log10(0.01/nrow(qtls))),]

getVarianceExplained <- function(genotypes, phenotypes, pheno.col = "mri42d_fat", marker = "UNC5048297"){
  ind           <- colnames(genotypes[marker,!is.na(genotypes[marker,])])

  genotype      <- as.factor(t(genotypes[marker,!is.na(genotypes[marker,])]))
  littersize    <- as.factor(phenotypes[ind, "WG2"])
  subfamily     <- as.factor(phenotypes[ind, "Vater"])
  season        <- as.factor(phenotypes[ind, "Season"])
  litternumber  <- as.factor(phenotypes[ind, "W.Label"])

  phenotype     <- phenotypes[ind, pheno.col]
  model <- lm(phenotype ~ subfamily + littersize + litternumber + season + genotype)
  tryCatch(res  <- anova(model), error = function(e){ res <<- NA })
  cat(names(model$coefficients),"\n")
  cat(model$coefficients,"\n")
  varExplained  <- res[, "Sum Sq"] / sum((phenotype - mean(phenotype, na.rm=TRUE))^2, na.rm=TRUE)
  names(varExplained) <- c("subfamily","littersize","litternumber", "season", "marker", "Left")
  return(round(varExplained * 100, digits=1))
}

topMarkerPerChromosome <- function(qtls, map, pheno.col="qtl42"){
  topmarkers <- NULL
  for(chr in unique(map[,"Chr"])){
    markers <- rownames(map[which(map[,"Chr"] == chr),])
    topmarkers <- c(topmarkers, rownames(qtls[markers,][which.max(qtls[markers,pheno.col]),])[1])
  }
  names(topmarkers) <- unique(map[,"Chr"])
  return(topmarkers)
}

mriGWAS_Cof <- function(genotypes, phenotypes, pheno.col = "42d", to = nrow(genotypes)){                                            # Use the chromosome 3 locus as cofactor
  pvalues <- NULL
  for(x in 1:to){
    ind           <- colnames(genotypes[x,!is.na(genotypes[x,])])
    genotype      <- as.factor(t(genotypes[x,!is.na(genotypes[x,])]))
    cof           <- as.factor(t(genotypes["UNC5048297", ind]))
    littersize    <- as.numeric(phenotypes[ind, "WG2"])
    phenotype     <- phenotypes[ind, paste0("mri",pheno.col,"_fat")] / phenotypes[ind, paste0("mri",pheno.col,"_lean")]
    tryCatch(res  <- anova(lm(phenotype ~ littersize + cof + genotype))[[5]][3], error = function(e){ res <<- NA })
    pvalues <- c(pvalues, res)
  }
  plot(-log10(pvalues),t='l')
  names(pvalues) <- rownames(genotypes)[1:to]
  return(-log10(pvalues))
}
qtl42cof   <- mriGWAS_Cof(genotypes, phenotypes, "42d")
plot(qtl42cof, t='h', col=chrcolors[map[,"Chr"]],main="Fat/Lean QTL profile Day 70", xlab="Marker", ylab="-log10(p-value)")

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))
setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
mlength      <- max(chrInfo[,"Length"])

plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="ASE & Dominant expression (maternal BFMI)", yaxt="n", xlab="Chromosome", ylab="Length (Mb)", xaxt="n")
abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

aa <- lapply(rownames(qtls), function(x){
  xloc <- match(as.character(map[x,"Chr"]), chromosomes); yloc <- as.numeric(map[x,"Mb_NCBI38"])
  cat(x, yloc, xloc,"\n")
  LOD <- qtls[x,1]
  points(x=xloc + (LOD/100), y=yloc , pch=15,cex=0.5)
})

cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1,lwd=3)
  cnt <<- cnt + 1
})

axis(1, chromosomes, at=c(1:nrow(chrInfo)), las=1)
axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
legend("topright", c("BFMI like", "B6N like"), fill=c("orange","gray"))


region1 <- "3:36481201:36854743"
regionO <- "3:34000000:44000000"
region2 <- "3:14841429:17563072"

library(biomaRt)
setwd("E:/Mouse/RNA/ArrayDesign/Atlas data")

mart      <- useMart("ensembl", "mmusculus_gene_ensembl")
allgenes  <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), filter=c("chromosomal_region"), values=regionO, mart = mart)


#### Plot chromosome 3, and looking into the second QTL

setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
qtls <- read.table("Analysis/qtls_mri70d_fatDlean.txt",   sep="\t")

ch3 <- rownames(map[which(map[rownames(qtls),1]==3),])
plot(qtls[ch3,],t='l',lwd=2,main="Fat / Lean")


qtls <- read.table("Analysis/qtls_mri70d_fatDlean_cof_UNC5048297.txt",   sep="\t")
qtls["UNC4768063",]

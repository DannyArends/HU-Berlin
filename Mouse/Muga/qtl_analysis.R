# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

getSeason <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES),".", fixed=TRUE),"[",2)))
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5] <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8] <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11] <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2] <- "Winter"
  return(ret)
}

setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                              # Normal A, H, B genotypes
genotypesPh <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("","AA","CC","TT","GG"))        # Phased by Beagle (only heterozygous)
genotypesGP <- read.table("Analysis/genotypesPhasedGP.txt",     sep="\t", check.names=FALSE, colClasses="character", na.strings=c(""))                            # Phased towards the grandparents

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                                                                    # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("Vater", "W.dat", "W.Label", "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
            
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)), phenos]                  # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F2 <- F2[-which(F2=="6661459")]                                                                               # This individual has no genotype data

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
birthmonth <- unlist(lapply(strsplit(as.character(phenotypes[,"W.dat"]),".", fixed=TRUE),"[",2))
phenotypes <- cbind(phenotypes, Birthmonth = birthmonth)                                                      # Add the birth month column to the matrix


onegenotype <- which(lapply(apply(genotypes[,F2], 1, table), length) == 1)                                    # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- unique(genotypes[-onegenotype, F2])                                                            # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypes), "markers\n")                                                                # == Left with 11677 markers

onegenotype <- which(lapply(apply(genotypesPh[,F2], 1, table), length) == 1)                                  # Markers with only one genotype cannot be used in QTL mapping
genotypesPh <- unique(genotypesPh[-onegenotype,F2])                                                           # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesPh), "markers\n")                                                              # == Left with 19468 markers

onegenotype <- which(lapply(apply(genotypesGP[,F2], 1, table), length) == 1)                                  # Markers with only one genotype cannot be used in QTL mapping
genotypesGP <- unique(genotypesGP[-onegenotype,F2])                                                           # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesGP), "markers\n")                                                              # == Left with 7585 markers

mriGWAS <- function(genotypes, phenotypes, pheno.col = "42d", to = nrow(genotypes)){
  pvalues <- NULL
  for(x in 1:to){
    ind            <- colnames(genotypes[x,!is.na(genotypes[x,])])                                                            # Which individuals have genotype data

    subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure (factor)
    littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter  (linear effect)
    litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter    (factor)
    season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born    (factor)
    genotype       <- as.factor(t(genotypes[x,!is.na(genotypes[x,])]))                                                        # The genotype under investigation  (factor)

    phenotype      <- phenotypes[ind, paste0("mri",pheno.col,"_fat")] / phenotypes[ind, paste0("mri",pheno.col,"_lean")]      # Response: Fat / Lean

    tryCatch(res <- anova(lm(phenotype ~ subfamily + littersize + litternumber + season + genotype + littersize:litternumber))[[5]], error = function(e){ res <<- rep(NA,5) })
    cat(x, round(-log10(res[-length(res)]),1),"\n")
    pvalues <- rbind(pvalues, res[-length(res)])
  }

  colnames(pvalues) <- c("subfamily", "l_size", "l_number", "season", "marker", "size:number")
  rownames(pvalues) <- rownames(genotypes)[1:to]
  return(round(-log10(pvalues), 3))
}

qtl42   <- mriGWAS(genotypes,   phenotypes, "42d") ; qtl56   <- mriGWAS(genotypes,   phenotypes, "56d") ; qtl70   <- mriGWAS(genotypes,   phenotypes, "70d")
qtlPH42 <- mriGWAS(genotypesPh, phenotypes, "42d") ; qtlPH56 <- mriGWAS(genotypesPh, phenotypes, "56d") ; qtlPH70 <- mriGWAS(genotypesPh, phenotypes, "70d")
qtlGP42 <- mriGWAS(genotypesGP, phenotypes, "42d") ; qtlGP56 <- mriGWAS(genotypesGP, phenotypes, "56d") ; qtlGP70 <- mriGWAS(genotypesGP, phenotypes, "70d")

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
write.table(qtl42,   "Analysis/qtls_fatDlean42.txt",   sep="\t")
write.table(qtl56,   "Analysis/qtls_fatDlean56.txt",   sep="\t")
write.table(qtl70,   "Analysis/qtls_fatDlean70.txt",   sep="\t")

write.table(qtlPH42, "Analysis/qtls_fatDlean_gwasPH.txt", sep="\t")
write.table(qtlPH56, "Analysis/qtls_fatDlean_gwasPH.txt", sep="\t")
write.table(qtlPH70, "Analysis/qtls_fatDlean_gwasPH.txt", sep="\t")

write.table(cbind(qtlGP42, qtlGP56, qtlGP70), "Analysis/qtls_fatDlean_gwasGP.txt", sep="\t")

qtls <- read.table("Analysis/qtls_fatDlean_gwas.txt", sep="\t", colClasses=c("character",rep("numeric",3)), header=TRUE)

image(x = 1:nrow(qtls), y=(1:3)-0.5, as.matrix(qtls), oldstyle=TRUE, breaks=c(0,3,6,12,30,100), col=c("white",gray.colors(4)[4:1]))
grid(3); box()

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
chrcolors <- rep(c("black","orange"),length(unique(map[,"Chr"])))
names(chrcolors) <- unique(map[,"Chr"])

png("Plots/QTL_FatdLean_Day42.png", width=1024, height=768)
  plot(qtls[,"qtl42"], t='h', col=chrcolors[map[,"Chr"]],main="Fat/Lean QTL profile Day 42", xlab="Marker", ylab="-log10(p-value)")
  abline(h = -log10(0.1/nrow(qtls)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtls)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtls)), col="green", lty=2)
dev.off()

png("Plots/QTL_FatdLean_Day56.png", width=1024, height=768)
  plot(qtls[,"qtl56"], t='h', col=chrcolors[map[,"Chr"]],main="Fat/Lean QTL profile Day 56", xlab="Marker", ylab="-log10(p-value)")
  abline(h = -log10(0.1/nrow(qtls)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtls)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtls)), col="green", lty=2)
dev.off()

png("Plots/QTL_FatdLean_Day70.png", width=1024, height=768)
  plot(qtls[,"qtl70"], t='h', col=chrcolors[map[,"Chr"]],main="Fat/Lean QTL profile Day 70", xlab="Marker", ylab="-log10(p-value)")
  abline(h = -log10(0.1/nrow(qtls)), col="orange", lty=2); abline(h = -log10(0.05/nrow(qtls)), col="gold", lty=2); abline(h = -log10(0.01/nrow(qtls)), col="green", lty=2)
dev.off()

map[which(qtls[,"qtl42"] > -log10(0.01/nrow(qtls))),]
map[which(qtls[,"qtl56"] > -log10(0.01/nrow(qtls))),]
map[which(qtls[,"qtl70"] > -log10(0.01/nrow(qtls))),]

getVarianceExplained <- function(genotypes, phenotypes, pheno.col = "42d", marker = "UNC5048297"){
  ind           <- colnames(genotypes[marker,!is.na(genotypes[marker,])])
  genotype      <- as.factor(t(genotypes[marker,!is.na(genotypes[marker,])]))
  littersize    <- as.factor(phenotypes[ind, "WG2"])
  subfamily     <- as.factor(phenotypes[ind, "Vater"])
  season        <- as.factor(phenotypes[ind, "Season"])
  litternumber  <- as.factor(phenotypes[ind, "W.Label"])
  phenotype     <- phenotypes[ind, paste0("mri",pheno.col,"_fat")] / phenotypes[ind, paste0("mri",pheno.col,"_lean")]
  tryCatch(res  <- anova(lm(phenotype ~ littersize + litternumber + subfamily + season + genotype + littersize:litternumber)), error = function(e){ res <<- NA })
  varExplained  <- res[, "Sum Sq"] / sum((phenotype - mean(phenotype, na.rm=TRUE))^2, na.rm=TRUE)
  names(varExplained) <- c("l_size","l_number","subfamily", "season", "marker", "Int", "Left")
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

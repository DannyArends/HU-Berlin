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

setwd("E:/Mouse/DNA/MegaMuga/")                                                                               # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)                              # Normal A, H, B genotypes
genotypesPh <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE)                  # Phased by Beagle, TODO: Leave out the homozygous individuals
genotypesGP <- read.table("Analysis/genotypesPhasedGP.txt", sep="\t", check.names=FALSE, na.strings="")       # Phased towards the grandparents

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("Vater", "W.dat", "W.Label", "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
            
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)), phenos]                  # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F2 <- F2[-which(F2=="6661459")]                                                                               # This individual has no genotype data

updown <- apply(phenotypes[F2,c("mri42d_fat", "mri56d_fat", "mri70d_fat")], 1, function(x){ sign(cor(1:3, as.numeric(x))) })

loosingFat <- rep(NA, length(rownames(phenotypes)))
names(loosingFat) <- rownames(phenotypes)
loosingFat[names(updown)] <- updown
phenotypes <- cbind(phenotypes, loosingFAT = loosingFat)

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))

genotypes   <- genotypes[,F2]
genotypesPh <- genotypesPh[,F2]
genotypesGP <- genotypesGP[,F2]

mriGWAS <- function(genotypes, phenotypes, pheno.col = "42d", to = nrow(genotypes)){                          # TODO: Add to the model: subfamily based on which F1 they come from
  pvalues <- NULL                                                                                             # TODO: ASK SEBASTIAAN: Add to the model: season, when were they born
  for(x in 1:to){                                                                                             # TODO: ASK SEBASTIAAN: Add to the model: litter number (1st litter versus second litter)
    ind            <- colnames(genotypes[x,!is.na(genotypes[x,])])
    genotype       <- as.factor(t(genotypes[x,!is.na(genotypes[x,])]))
    littersize     <- as.numeric(phenotypes[ind, "WG2"])
    subfamily      <- as.factor(phenotypes[ind, "Vater"])
    season         <- as.factor(phenotypes[ind, "Season"])
    litternumber   <- as.factor(phenotypes[ind, "W.Label"])
    #phenotype     <- phenotypes[ind, paste0("mri",pheno.col,"_fat")] / phenotypes[ind, paste0("mri",pheno.col,"_lean")]
    phenotype     <- phenotypes[ind, pheno.col]
    tryCatch(res  <- anova(lm(phenotype ~ littersize + litternumber + subfamily + season + genotype))[[5]], error = function(e){ res <<- rep(NA,5) })
    cat(x, res,"\n")
    pvalues <- c(pvalues, res[5])
  }
  plot(-log10(pvalues),t='l')
  names(pvalues) <- rownames(genotypes)[1:to]
  return(-log10(pvalues))
}

qtlFAT  <- mriGWAS(genotypes,   phenotypes, "loosingFAT")

qtl42   <- mriGWAS(genotypes,   phenotypes, "42d") ; qtl56   <- mriGWAS(genotypes,   phenotypes, "56d") ; qtl70   <- mriGWAS(genotypes,   phenotypes, "70d")
qtlPH42 <- mriGWAS(genotypesPh, phenotypes, "42d") ; qtlPH56 <- mriGWAS(genotypesPh, phenotypes, "56d") ; qtlPH70 <- mriGWAS(genotypesPh, phenotypes, "70d")
qtlGP42 <- mriGWAS(genotypesGP, phenotypes, "42d") ; qtlGP56 <- mriGWAS(genotypesGP, phenotypes, "56d") ; qtlGP70 <- mriGWAS(genotypesGP, phenotypes, "70d")

setwd("E:/Mouse/ClassicalPhenotypes/AIL")

write.table(cbind(qtl42,   qtl56,   qtl70),   "Analysis/qtls_fatDlean_gwas.txt",   sep="\t")
write.table(cbind(qtlPH42, qtlPH56, qtlPH70), "Analysis/qtls_fatDlean_gwasPH.txt", sep="\t")
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
  littersize    <- as.numeric(phenotypes[ind, "WG2"])
  phenotype     <- phenotypes[ind, paste0("mri",pheno.col,"_fat")] / phenotypes[ind, paste0("mri",pheno.col,"_lean")]
  tryCatch(res  <- anova(lm(phenotype ~ littersize + genotype)), error = function(e){ res <<- NA })
  varExplained  <- res[2, "Sum Sq"] / sum((phenotype - mean(phenotype, na.rm=TRUE))^2, na.rm=TRUE)
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

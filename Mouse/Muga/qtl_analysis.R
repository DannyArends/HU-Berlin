# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/DNA/MegaMuga/")                                                                               # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)                              # Normal A, H, B genotypes
genotypesPh <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE)                  # Phased by Beagle
genotypesGP <- read.table("Analysis/genotypesPhasedGP.txt", sep="\t", check.names=FALSE, na.strings="")       # Phased towards the grandparents

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),phenos]                   # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F2 <- F2[-which(F2=="6661459")]                                                                               # This individual has no genotype data

genotypes   <- genotypes[,F2]
genotypesPh <- genotypesPh[,F2]
genotypesGP <- genotypesGP[,F2]

mriGWAS <- function(genotypes, phenotypes, pheno.col = "42d", to = nrow(genotypes)){
  pvalues <- NULL
  for(x in 1:to){
    ind           <- colnames(genotypes[x,!is.na(genotypes[x,])])
    genotype      <- as.factor(t(genotypes[x,!is.na(genotypes[x,])]))
    littersize    <- as.numeric(phenotypes[ind, "WG2"])
    phenotype     <- phenotypes[ind, paste0("mri",pheno.col,"_fat")] / phenotypes[ind, paste0("mri",pheno.col,"_lean")]
    tryCatch(res  <- anova(lm(phenotype ~  littersize + genotype))[[5]][2], error = function(e){ res <<- NA })
    pvalues <- c(pvalues, res)
  }
  plot(-log10(pvalues),t='l')
  names(pvalues) <- rownames(genotypes)[1:to]
  return(-log10(pvalues))
}

qtl42   <- mriGWAS(genotypes,   phenotypes, "42d") ; qtl56   <- mriGWAS(genotypes,   phenotypes, "56d") ; qtl70   <- mriGWAS(genotypes,   phenotypes, "70d")
qtlPH42 <- mriGWAS(genotypesPh, phenotypes, "42d") ; qtlPH56 <- mriGWAS(genotypesPh, phenotypes, "56d") ; qtlPH70 <- mriGWAS(genotypesPh, phenotypes, "70d")
qtlGP42 <- mriGWAS(genotypesGP, phenotypes, "42d") ; qtlGP56 <- mriGWAS(genotypesGP, phenotypes, "56d") ; qtlGP70 <- mriGWAS(genotypesGP, phenotypes, "70d")

setwd("E:/Mouse/ClassicalPhenotypes/AIL")

write.table(cbind(qtl42,   qtl56,   qtl70), "Analysis/qtls_fatDlean_gwas.txt", sep="\t")
write.table(cbind(qtlPH42, qtlPH56, qtlPH70), "Analysis/qtls_fatDlean_gwasPH.txt", sep="\t")
write.table(cbind(qtlGP42, qtlGP56, qtlGP70), "Analysis/qtls_fatDlean_gwasGP.txt", sep="\t")

qtls <- read.table("Analysis/qtls_fatDlean_gwas.txt", sep="\t", colClasses=c("character",rep("numeric",3)), header=TRUE,row.names=FALSE)

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

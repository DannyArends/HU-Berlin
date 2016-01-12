
source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

# Load the genotype call data
setwd("E:/Mouse/DNA/MegaMuga/")
phased.vcf <- read.table(gzfile(paste0("Analysis/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("Analysis/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header

samples <- colnames(phased.vcf)[-c(1:9)]

phased.vcf[, samples] <- apply(phased.vcf[, samples],2,function(x){unlist(lapply(strsplit(x,":"),"[",1))})                   # Only the genotypes
rownames(phased.vcf) <- phased.vcf[,"ID"]
phased.vcf <- phased.vcf[,-which(colnames(phased.vcf) %in% c("ID", "QUAL","FILTER","INFO","FORMAT"))]                              # Remove unneeded columns
phased.vcf[1:10,1:10]

# Load in the phenotypes and create the pedigree file
phenotypes <- read.table("Phenotypes/allPhenotypes.txt",sep="\t",header=TRUE,na.strings=c(0,"-","NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.AHBp)),]                # We do not have genotypes for all individuals

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals

# Change the genotype coding
phased.AHBp <- phased.vcf                                                                         # Copy
phased.AHBp[, samples] <- apply(phased.AHBp[, samples], 2, fromVCFcode.AHBp)                      # Change coding to A H0 H1 B


# QTL Analysis of generation 28
# Filter the genotypes for markers that are seggregating
phased.AHBp <- phased.AHBp[-which(lapply(apply(phased.AHBp[, F2], 1, table),length) == 1),]

# Create some H0 versus H1 statistics
statistics <- apply(phased.AHBp[, F2], 2, table)
plot(statistics[3,], statistics[4,], xlab="H0 - A from father", ylab="H1 - B from father", main="H0 versus H1 in Generation 28")

# Now we want to do the QTL mapping for a single phenotype
covariates <- phenotypes[F2, c("Eltern_ID", "WG2", "W.Label", "Season")]
phenotype  <- phenotypes[F2, "mri70d_fat"]

cnt <- 1
pvalues <- t(apply(phased.AHBp[,F2], 1, function(x){
  mylm <- lm(phenotype ~ covariates[,"Eltern_ID"] + covariates[,"WG2"] + covariates[,"W.Label"] + covariates[,"Season"] + as.factor(as.character(x)))
  cnt <<- cnt + 1
  unlist(anova(mylm)[[5]])
}))
plot(-log10(p.adjust(pvalues[,5],"BH")), col = as.numeric(as.factor(phased.AHBp[,"CHROM"])),t ='h')


ratios <- apply(phased.AHBp[,F2],1,function(x){
  return(sum(x == "H0") / (sum(x == "H0")+sum(x == "H1")))
})




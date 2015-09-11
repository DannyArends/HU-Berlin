
toVCFcode <- function(x){
  r <- rep(NA,length(x));  r[x == "A"] <- "0/0"; r[x == "H"] <- "0/1";  r[x == "B"] <- "1/1";  r[is.na(x)] <- "./.";  return(r)
}

fromVCFcode.AHB <- function(x){
  r <- rep(NA,length(x));  r[x == "0|0"] <- "A";  r[x == "0|1"] <- "H"; r[x == "1|0"] <- "H";  r[x == "1|1"] <- "B";  r[x == ".|."] <- NA;  return(r)
}

fromVCFcode.AHBp <- function(x){
  r <- rep(NA,length(x));  r[x == "0|0"] <- "A";  r[x == "0|1"] <- "H0"; r[x == "1|0"] <- "H1";  r[x == "1|1"] <- "B";  r[x == ".|."] <- NA;  return(r)
}

# Load the genotype call data
setwd("E:/Mouse/DNA/MegaMuga/")
locusxdnaheader     <- unlist(strsplit(readLines("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", n=16)[16],","))
locusxdna           <- read.csv("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", header=FALSE, skip = 22)
colnames(locusxdna) <- c("Label","plateWell","Date","oligoPoolId","bundleId", "status", "Type", "Nas", locusxdnaheader[4:length(locusxdnaheader)])
locusxdna           <- t(locusxdna[which(locusxdna[,"Type"] == "calls"),])
colnames(locusxdna) <- locusxdna[1,]
locusxdna <- locusxdna[-c(1:8), ]                                                     # Remove some unneeded headers
locusxdna[which(locusxdna == "U")] <- NA                                              # Use the correct NA values

# Remove markers and individuals without any data and repeated measurements for the same marker
locusxdna <- locusxdna[-which(apply(locusxdna,1,function(x){return(sum(is.na(x)))}) == ncol(locusxdna)),]
locusxdna <- locusxdna[,-which(apply(locusxdna,2,function(x){return(sum(is.na(x)))}) == nrow(locusxdna))]

cat("Loaded:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 479 individuals, 70651 markers
write.table(locusxdna, "Analysis/AllGenotypes.txt", sep="\t", quote=FALSE)

samples <- colnames(locusxdna)
markers <- rownames(locusxdna)

#Load annotation from JAX
load("MM_snps.Rdata")
MM_snps[,"SNP_ID"] <- gsub(".", "-", MM_snps[,"SNP_ID"], fixed = TRUE)                # Fix the . to - error in the MM_SNP
rownames(MM_snps) <- MM_snps[,"SNP_ID"]                                               # Use them as rownames

marker.annot <- matrix(NA, length(markers), 4, dimnames=list(markers,c("Chr", "Pos", "snp_A", "snp_B")))
marker.annot[markers, "Chr"] <- MM_snps[markers,"Chr"]
marker.annot[markers, "Pos"] <- MM_snps[markers,"Mb_NCBI38"] * 1000000

for(x in markers){
  instr  <- regexpr("\\[.+\\]", MM_snps[x,"Sequence"], perl=TRUE)
  marker.annot[x,c("snp_A","snp_B")] <- unlist(strsplit(substr(MM_snps[x,"Sequence"], (instr[1]+1), (instr[1]+attr(instr,"match.length")-2)),"/"))
}
write.table(marker.annot, "Analysis/markerAnnotation.txt", sep="\t", quote=FALSE)

# Remove the duplicated markers
duplicates <- rownames(marker.annot)[which(duplicated(apply(marker.annot[,c("Chr","Pos")],1,paste0,collapse="")))]
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% duplicates),]
locusxdna <- locusxdna[-which(rownames(locusxdna) %in% duplicates),]

cat("After duplicate removal:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 479 individuals, 70651 markers

# Remove markers that do not have a known location
missingloc <- rownames(marker.annot)[which(is.na(marker.annot[,"Pos"]))]
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% missingloc),]
locusxdna <- locusxdna[-which(rownames(locusxdna) %in% missingloc),]

cat("After duplicate and unplaced marker removal:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 479 individuals, 70651 markers

## Order the markers since the VCF needs to have continuous markers per chromosome
marker.annot <- marker.annot[order(marker.annot[,"Chr"], as.numeric(marker.annot[,"Pos"])),]
markers <- rownames(marker.annot)
locusxdna <- locusxdna[markers, ]

# Write genotypes as VCF for Beagle4
if(!file.exists("Analysis/genotypes.vcf")){
  cat("##fileformat=VCFv4.2\n", file="Analysis/genotypes.vcf")
  cat("##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n", file="Analysis/genotypes.vcf", append = TRUE)
  cat(paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste(samples, collapse = "\t"), "\n"), file="Analysis/genotypes.vcf", append = TRUE)
  for(x in markers){
    gts <- paste(toVCFcode(locusxdna[x, samples]), collapse = "\t")
    cat(marker.annot[x,"Chr"], "\t", marker.annot[x,"Pos"], "\t", x, "\t", marker.annot[x,"snp_A"], "\t", marker.annot[x,"snp_B"], "\t.\tPASS\t.\tGT\t", gts, "\n", sep="", file="Analysis/genotypes.vcf", append = TRUE)
  }
}

# Load in the phenotypes and create the pedigree file
phenotypes <- read.table("Phenotypes/allPhenotypes.txt",sep="\t",header=TRUE,na.strings=c(0,"-","NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix

if(!file.exists("Analysis/population.ped")){
  cat("", file="Analysis/population.ped") # Empty ped file holding the pedigree
  for(x in samples){
    dataline <- c("1", x, phenotypes[x,"Vater"], phenotypes[x,"Mutter"])
    dataline[is.na(dataline)] <- "0"
    cat(paste(dataline, collapse=" "), "\n", sep = "", file="Analysis/population.ped", append=TRUE) # Empty ped file holding the pedigree
  }
}

if(!file.exists("Analysis/phased.vcf.gz")){ # Do beagle phasing
  system(paste0("java -Xmx1000m -jar beagle.v.4.jar gt=Analysis/genotypes.vcf ped=Analysis/population.ped out=Analysis/phased"))
}

# After phasing load the created VCF file with genotypes
phased.vcf <- read.table(gzfile(paste0("Analysis/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("Analysis/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header
phased.vcf[, samples] <- apply(phased.vcf[, samples],2,function(x){unlist(lapply(strsplit(x,":"),"[",1))})                   # Only the genotypes
rownames(phased.vcf) <- phased.vcf[,"ID"]
phased.vcf <- phased.vcf[,-which(colnames(phased.vcf) %in% c("ID", "QUAL","FILTER","INFO","FORMAT"))]                              # Remove unneeded columns
phased.vcf[1:10,1:10]

# Change the genotype coding
phased.AHBp <- phased.vcf                                                                         # Copy
phased.AHBp[, samples] <- apply(phased.AHBp[, samples], 2, fromVCFcode.AHBp)                      # Change coding to A H0 H1 B

# Order chromosomes in the normal way
chromosomes  <- c(1:19, "X", "Y", "M")
phased.AHBpN <- NULL
for(chr in chromosomes){
  phased.AHBpN <- rbind(phased.AHBpN, phased.AHBp[which(phased.AHBp[,"CHROM"] == chr),])
}
phased.AHBp <- phased.AHBpN

# Load in the phenotypes
phenotypes <- read.table("Phenotypes/allPhenotypes.txt", sep="\t", header=TRUE, na.strings=c(0, "-", "NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                                         # Add the season column to the matrix
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.AHBp)),]                                               # We do not have genotypes for all individuals

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)] ; F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]     # Get the names for each generations

# parent of origin effects

for(m in 1:nrow(phased.AHBp)){
  pTrans <- matrix(0,2,2,dimnames=list(c("A","B"),c("A","B")))
  for(i in F2){
    pG <- phased.AHBp[m, as.character(phenotypes[i,"Vater"])]
    iG <- phased.AHBp[m, i]
    if(pG == "A" || pG == "B")   pTrans[pG,pG] <- pTrans[pG,pG] + 1
    if(pG == "H0" || pG == "H1"){
      if(iG == "A") pTrans["A","B"] <- pTrans["A","B"] + 1
      if(iG == "B") pTrans["B","A"] <- pTrans["B","A"] + 1
    }
    if(iG == "H0") 
    if(iG == "H1")
  }
}


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




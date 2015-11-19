

source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
source("D:/Github/HU-Berlin/Mouse/Muga/vcfTools.R")

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
if(!file.exists("Analysis/AllGenotypes.txt")) write.table(locusxdna, "Analysis/AllGenotypes.txt", sep="\t", quote=FALSE)

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


F1m <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27 & phenotypes[, "sex"] == 'm')]
F1f <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27 & phenotypes[, "sex"] == 'f')]

counts <- matrix(c(rep(0,8),rep(NA,4)), nrow(phased.AHBp), 12, dimnames = list(rownames(phased.AHBp), c("P12","P21", "M12", "M21", "C12", "C21", "nFat", "nMat", "Xp", "Xm", "PO", "Xh")),byrow=TRUE)

updateCounts <- function(m, type, pG, iG){
  if(pG == "H0" || pG == "H1"){
    if(iG == "A")  counts[m, paste0(type,"21")] <<- counts[m, paste0(type,"21")] + 1
    if(iG == "B")  counts[m, paste0(type,"12")] <<- counts[m, paste0(type,"12")] + 1
    if(iG == "H0") counts[m, paste0(type,"21")] <<- counts[m, paste0(type,"21")] + 1
    if(iG == "H1") counts[m, paste0(type,"12")] <<- counts[m, paste0(type,"12")] + 1
  }
}

for(m in 1:nrow(phased.AHBp)){
  for(i in F2){
    vG <- phased.AHBp[m, as.character(phenotypes[i, "Vater"])]
    mG <- phased.AHBp[m, as.character(phenotypes[i, "Mutter"])]
    iG <- phased.AHBp[m, i]
    updateCounts(m, "P", vG, iG)
    updateCounts(m, "M", mG, iG)
    if(vG == "H0" || vG == "H1" ){
      if(mG == "H0" || mG == "H1" ){
        if(iG == "H0") counts[m, "C12"] <- counts[m, "C12"] + 1
        if(iG == "H1") counts[m, "C21"] <- counts[m, "C21"] + 1
      }
    }
  }
  counts[m, "nFat"] <-  length(which(grepl("H", phased.AHBp[m,F1m])))
  counts[m, "nMat"] <-  length(which(grepl("H", phased.AHBp[m,F1f])))

  counts[m, "Xp"] <-  ((counts[m, "P12"] - counts[m, "P21"])^2) / (counts[m, "P12"] + counts[m, "P21"])
  counts[m, "Xm"] <-  ((counts[m, "M12"] - counts[m, "M21"])^2) / (counts[m, "M12"] + counts[m, "M21"])
  counts[m, "PO"] <-  log(counts[m, "P12"]/counts[m, "P21"]) - log(counts[m, "M12"]/counts[m, "M21"]) / sqrt((1/counts[m, "P12"]) + (1/counts[m, "P21"]) + (1/counts[m, "M12"]) + (1/counts[m, "M21"]))
  counts[m, "Xh"] <-  ((counts[m, "C12"] - counts[m, "C21"])^2) / (counts[m, "C12"] + counts[m, "C21"])
  cat(m,"/",nrow(phased.AHBp),"\n")
}

write.table(counts, file="Analysis/TransmissionBias.txt", sep="\t")

counts <- read.table("Analysis/TransmissionBias.txt", sep="\t", row.names=1, header=TRUE)

enoughSamples <- unlist(apply(counts[,c("P12", "P21", "M12", "M21")] ,1,function(x){
  all(x > 5)
}))

counts <- counts[enoughSamples,]
map <- marker.annot[rownames(counts), ]
map.autosomes <- map[which(map[,"Chr"] %in% 1:19),]

counts <- counts[rownames(map.autosomes),]

Pp  <- -log10(pchisq(counts[, "Xp"], 1, lower.tail=FALSE))
names(Pp) <- rownames(map.autosomes)
Pm  <- -log10(pchisq(counts[, "Xm"], 1, lower.tail=FALSE))
names(Pm) <- rownames(map.autosomes)
PoO <-  -log10(pchisq(counts[, "PO"], 1, lower.tail=FALSE))
Ph  <-  -log10(pchisq(counts[, "Xh"], 1, lower.tail=FALSE))

bfmi <- phased.AHBp[rownames(counts[which(Pp > 10),]), which(colnames(phased.AHBp) == "BFMI860-12 (V2)")]
patTransmission <- cbind(counts[rownames(counts[which(Pp > 10),]),],bfmi)

transmitted <- NULL
for(x in 1:nrow(patTransmission)){
  if(patTransmission[x, "P12"] > patTransmission[x, "P21"]){ # prefere to transmit a B`
  transmitted <- c(transmitted, "B")
  }else{
  transmitted <- c(transmitted, "A")
  }
}

write.table(map.autosomes[rownames(counts[which(Pp > 10),]),],"Analysis/SignPp.txt",sep="\t")
write.table(cbind(map.autosomes[rownames(counts[which(Pp > 10),]),], patTransmission,transmitted, T = transmitted == bfmi),"Analysis/SignPpCounts.txt",sep="\t")

cat(length(which(bfmi == transmitted)), "BFMI out of", length(transmitted),"\n")

bfmi <- phased.AHBp[rownames(counts[which(Pm > 10),]), which(colnames(phased.AHBp) == "BFMI860-12 (V2)")]
matTransmission <- cbind(counts[rownames(counts[which(Pm > 10),]),],bfmi)

transmitted <- NULL
for(x in 1:nrow(matTransmission)){
  if(matTransmission[x, "P12"] > matTransmission[x, "P21"]){ # prefere to transmit a B`
  transmitted <- c(transmitted, "B")
  }else{
  transmitted <- c(transmitted, "A")
  }
}

write.table(map.autosomes[rownames(counts[which(Pm > 10),]),],"Analysis/SignPm.txt",sep="\t")
write.table(cbind(map.autosomes[rownames(counts[which(Pm > 10),]),], matTransmission,transmitted, T = transmitted == bfmi),"Analysis/SignPmCounts.txt",sep="\t")

cat(length(which(bfmi == transmitted)), "BFMI out of", length(transmitted),"\n")

chr.lengths <- NULL
lsum <- 0
for(chr in unique(map.autosomes[,"Chr"])){
  onChr <- map.autosomes[map.autosomes[,"Chr"] == chr,]
  l <- max(as.numeric(onChr[,"Pos"]))
  chr.lengths <- rbind(chr.lengths, c(l, lsum, lsum + l))
  lsum <- lsum + l
}
rownames(chr.lengths) <- unique(map.autosomes[,"Chr"])

chr.cols <- 1+as.numeric(map.autosomes[,"Chr"]) %%2
chr.cols2 <- 1+(as.numeric(map.autosomes[,"Chr"]) %%2 == 0)

ymax <- max(c(Pp,Pm),na.rm=TRUE)
i <- 1
plot(c(0,max(chr.lengths[,3])), c(-ymax, ymax), t = 'n', ylab="-log10(p-value)", xlab="Position (Mb)")
for(chr in unique(map.autosomes[,"Chr"])){
  onChr <- rownames(map.autosomes[map.autosomes[,"Chr"] == chr,])
  points(x=as.numeric(map.autosomes[onChr,"Pos"]) + chr.lengths[chr, 2], y = Pp[onChr], t = 'h', col=c("Blue", "Gray")[i])
  points(x=as.numeric(map.autosomes[onChr,"Pos"]) + chr.lengths[chr, 2], y = -Pm[onChr], t = 'h', col=c("Gray", "Pink")[i])
  if(i == 2){ i <- 1; }else{ i <- 2; }
}
abline(h=-log10(0.05 / (4*nrow(counts))), col="orange")
abline(h=-log10(0.01 / (4*nrow(counts))), col="green")
abline(h=log10(0.05 / (4*nrow(counts))), col="orange")
abline(h=log10(0.01 / (4*nrow(counts))), col="green")



plot(abs(Pp-Pm), t='h', col=c("black", "Gray")[chr.cols2])
abline(h=-log10(0.05 / (4*nrow(counts))), col="orange")
abline(h=-log10(0.01 / (4*nrow(counts))), col="green")

plot(PoO, t='h', col=c("orange", "green")[chr.cols])

## Create the chromosome plot
op <- par(mfrow=c(1,1))
ymax <- max(as.numeric(map.autosomes[,"Pos"]))

plot(c(1,19), c(0, ymax), t = 'n', xlab="Chromosome", ylab="Position (Mbp)", xaxt='n', yaxt='n', main="Transmission bias from heterozygous parents")
axis(1, at=1:19, paste0("Chr ", 1:19), las=2, cex.axis=0.7)
axis(2, at=seq(0,ymax, 20000000), seq(0, ymax, 20000000) / 1000000, las=2, cex.axis=0.9)
for(x in 1:19){
  onChr <- rownames(map.autosomes[which(map.autosomes[,"Chr"] == as.character(x)),])
  colz <- as.numeric(-log10(pchisq(counts[onChr, "Xp"], 1, lower.tail=FALSE)) > 7)
  colz <- colz + as.numeric(-log10(pchisq(counts[onChr, "Xp"], 1, lower.tail=FALSE)) > 10) + 1

  colfunc <- c("white", "lightblue", "darkblue")
  points(rep(x-0.15,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col=colfunc[colz], cex=1.8)

  colz <- as.numeric(-log10(pchisq(counts[onChr, "Xm"], 1, lower.tail=FALSE)) > 7)
  colz <- colz + as.numeric(-log10(pchisq(counts[onChr, "Xm"], 1, lower.tail=FALSE)) > 10) + 1
  colfunc <- c("white", "pink", "red")
  points(rep(x+0.15,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col=colfunc[colz], cex=1.8)

  onC <- which(marker.annot[,"Chr"] == x)
  points(rep(x,length(onC)), as.numeric(marker.annot[onC,"Pos"]), pch="-",cex=1.5, col="gray")
  points(rep(x,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col='black',cex=1.5)
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




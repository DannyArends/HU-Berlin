#
# Analyse ATB in megaMuga data from multiple generations after beagle phasing
#

source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

# Load the genotype call data
setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
phased.vcf <- read.table(gzfile(paste0("Analysis/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("Analysis/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header

samples <- colnames(phased.vcf)[-c(1:9)]

phased.vcf[, samples] <- apply(phased.vcf[, samples],2,function(x){unlist(lapply(strsplit(x,":"),"[",1))})                   # Only the genotypes
rownames(phased.vcf) <- phased.vcf[,"ID"]
phased.vcf <- phased.vcf[,-which(colnames(phased.vcf) %in% c("ID", "QUAL","FILTER","INFO","FORMAT"))]                              # Remove unneeded columns
phased.vcf[1:10,1:10]

# Change the genotype coding
phased.AHBp <- phased.vcf                                                                         # Copy
phased.AHBp[, samples] <- apply(phased.AHBp[, samples], 2, fromVCFcode.AHBp)                      # Change coding to A H0 H1 B

# Change the genotype coding
phased.geno <- phased.vcf                                                                         # Copy
phased.geno[, samples] <- apply(phased.geno[, samples], 2, fromVCFcode.geno)                      # Change coding to AA AB BA BB

# Order chromosomes in the normal way
chromosomes  <- c(1:19, "X", "Y", "M")
phased.AHBpN <- NULL
phased.genoN <- NULL
for(chr in chromosomes){
  phased.AHBpN <- rbind(phased.AHBpN, phased.AHBp[which(phased.AHBp[,"CHROM"] == chr),])
  phased.genoN <- rbind(phased.genoN, phased.geno[which(phased.geno[,"CHROM"] == chr),])
}
phased.AHBp <- phased.AHBpN
phased.geno <- phased.genoN

#Load marker annotation
marker.annot <- read.table("Analysis/markerAnnotation.txt", sep="\t")
# Remove the duplicated markers
duplicates <- rownames(marker.annot)[which(duplicated(apply(marker.annot[,c("Chr","Pos")],1,paste0,collapse="")))]
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% duplicates),]

# Remove markers that do not have a known location
missingloc <- rownames(marker.annot)[which(is.na(marker.annot[,"Pos"]))]
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% missingloc),]

#order the chromosomes
marker.annot <- marker.annot[order(marker.annot[,"Chr"], as.numeric(marker.annot[,"Pos"])),]
marker.annot <- marker.annot[rownames(phased.geno), ]

# Load in the phenotypes
phenotypes <- read.table("Phenotypes/allPhenotypes.txt", sep="\t", header=TRUE, na.strings=c(0, "-", "NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                          # Add the season column to the matrix
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.AHBp)),]                # We do not have genotypes for all individuals

calculateATB <- function(phased.AHBp, phenotypes, generation = 28){
  individuals <- rownames(phenotypes)[which(phenotypes[, "Gen."] == generation)]
  if(!file.exists(paste0("Analysis/TransmissionBias_",generation,".txt"))){
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
      for(i in individuals){
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
      F1m <- rownames(phenotypes)[which(phenotypes[, "Gen."] == (generation-1) & phenotypes[, "sex"] == 'm')]
      F1f <- rownames(phenotypes)[which(phenotypes[, "Gen."] == (generation-1) & phenotypes[, "sex"] == 'f')]

      counts[m, "nFat"] <-  length(which(grepl("H", phased.AHBp[m,F1m])))
      counts[m, "nMat"] <-  length(which(grepl("H", phased.AHBp[m,F1f])))

      counts[m, "Xp"] <-  ((counts[m, "P12"] - counts[m, "P21"])^2) / (counts[m, "P12"] + counts[m, "P21"])
      counts[m, "Xm"] <-  ((counts[m, "M12"] - counts[m, "M21"])^2) / (counts[m, "M12"] + counts[m, "M21"])
      counts[m, "PO"] <-  log10(counts[m, "P12"]/counts[m, "P21"]) - log10(counts[m, "M12"]/counts[m, "M21"]) / sqrt((1/counts[m, "P12"]) + (1/counts[m, "P21"]) + (1/counts[m, "M12"]) + (1/counts[m, "M21"]))
      counts[m, "Xh"] <-  ((counts[m, "C12"] - counts[m, "C21"])^2) / (counts[m, "C12"] + counts[m, "C21"])
      if(m %% 100 == 0) cat(m,"/",nrow(phased.AHBp),"\n")
    }
    write.table(counts, file=paste0("Analysis/TransmissionBias_",generation,".txt"), sep="\t")
  }else{
    cat("Loading results from disk:", paste0("Analysis/TransmissionBias_",generation,".txt"), "\n")
    counts <- read.table(paste0("Analysis/TransmissionBias_", generation, ".txt"), sep="\t", row.names=1, header=TRUE)
  }
  return(counts)
}

getSignificant <- function(phased.AHBp, phenotypes, generation = 28){
  counts <- NULL
  if(!file.exists(paste0("Analysis/TransmissionBias_annotated_",generation,".txt"))){
    counts <- calculateATB(phased.AHBp, phenotypes, generation)
    enoughSamples <- unlist(apply(counts[,c("nFat", "nMat")] , 1, function(x){ any(x > 10) }))
    counts <- counts[enoughSamples,]
    map <- marker.annot[rownames(counts), ]
    map.autosomes <- map[which(map[,"Chr"] %in% 1:19),]

    counts  <- counts[rownames(map.autosomes),]
    pPat    <- pchisq(counts[, "Xp"], 1, lower.tail=FALSE)          # Allele bias Father
    #lodPat  <- -log10(pchisq(counts[, "Xp"], 1, lower.tail=FALSE))  # Allele bias Father
    pMat    <- pchisq(counts[, "Xm"], 1, lower.tail=FALSE)          # Allele bias Mother
    #lodMat  <- -log10(pchisq(counts[, "Xm"], 1, lower.tail=FALSE))  # Allele bias Mother
    pPoO    <- 2*pnorm(-abs(counts[, "PO"]))                        # Ratio Pat/Mat
    #lodPoO  <- -log10(2*pnorm(-abs(counts[, "PO"])))                # Ratio Pat/Mat
    pPh     <- pchisq(counts[, "Xh"], 1, lower.tail=FALSE)          # Overrepresentation of H0 versus H1
    #lodPh   <- -log10(pchisq(counts[, "Xh"], 1, lower.tail=FALSE))  # Overrepresentation of H0 versus H1
    bfmi    <- phased.AHBp[rownames(counts), which(colnames(phased.AHBp) == "BFMI860-12 (V2)")]
    counts  <- cbind(counts, pPat, pMat, pPoO, pPh, bfmi, patPref = NA, matPref = NA)

    for(x in 1:nrow(counts)){
      if(counts[x, "P12"] > counts[x, "P21"]){ counts[x,"patPref"] <- "B"; }else{ counts[x,"patPref"] <- "A"; }
      if(counts[x, "M12"] > counts[x, "M21"]){ counts[x,"matPref"] <- "B"; }else{ counts[x,"matPref"] <- "A"; }
    }
    counts <- cbind(map.autosomes[rownames(counts),],counts)
    write.table(counts, file=paste0("Analysis/TransmissionBias_annotated_",generation,".txt"), sep="\t")
  }else{
    counts <- read.table(paste0("Analysis/TransmissionBias_annotated_",generation,".txt"))
  }
  return(counts)
}

counts27 <- getSignificant(phased.AHBp, phenotypes, generation = 27)
counts28 <- getSignificant(phased.AHBp, phenotypes, generation = 28)

sPat27 <- counts27[which(-log10(counts27[,"pPat"]) > -log10(0.05/(70000*4))),]
if(!file.exists("Analysis/TransmissionBias_Pat_0.05_27.txt")) {
  write.table(sPat27, file="Analysis/TransmissionBias_Pat_0.05_27.txt", sep="\t")
}
sPat28 <- counts28[which(-log10(counts28[,"pPat"]) > -log10(0.01/(70000*4))),]
if(!file.exists("Analysis/TransmissionBias_Pat_0.01_28.txt")) {
  write.table(sPat28, file="Analysis/TransmissionBias_Pat_0.01_28.txt", sep="\t")
}
all(rownames(sPat27) %in% rownames(sPat28))

sMat27 <- counts27[which(-log10(counts27[,"pMat"]) > -log10(0.05/(70000*4))),]
if(!file.exists("Analysis/TransmissionBias_Mat_0.05_27.txt")) {
  write.table(sMat27, file="Analysis/TransmissionBias_Mat_0.05_27.txt", sep="\t")
}
sMat28 <- counts28[which(-log10(counts28[,"pMat"]) > -log10(0.01/(70000*4))),]
if(!file.exists("Analysis/TransmissionBias_Mat_0.01_28.txt")) {
  write.table(sMat28, file="Analysis/TransmissionBias_Mat_0.01_28.txt", sep="\t")
}
all(rownames(sMat27) %in% rownames(sMat28))

create.plot <- function(counts, HWEdata = NULL){
  map <- marker.annot[rownames(counts), ]
  map.autosomes <- map[which(map[,"Chr"] %in% 1:19),]

  chr.lengths <- NULL
  lsum <- 0
  for(chr in unique(map.autosomes[,"Chr"])){
    onChr <- map.autosomes[map.autosomes[,"Chr"] == chr,]
    l <- max(as.numeric(onChr[,"Pos"]))
    chr.lengths <- rbind(chr.lengths, c(l, lsum, lsum + l))
    lsum <- lsum + l
  }
  rownames(chr.lengths) <- unique(map.autosomes[,"Chr"])

  map.autosomes <- cbind(map.autosomes, cumPos = NA)

  for(chr in unique(map.autosomes[,"Chr"])){
    onChr <- rownames(map.autosomes[map.autosomes[,"Chr"] == chr,])
    map.autosomes[onChr,"cumPos"] <- map.autosomes[onChr,"Pos"] + chr.lengths[chr,2]
  }

  N <- 4*70000
  op <- par(mfrow = c(4 + as.numeric(!is.null(HWEdata)), 1), mar = c(3,4,1.5,1))

  chr.cols <- c("black","gray")[1+as.numeric(map.autosomes[,"Chr"]) %%2]

  plot(x= map.autosomes[,"cumPos"],  y = -log10(counts[,"pPat"]), col = chr.cols, pch=19,xaxt='n', ylab="-log10(Ppat)", main = "PAT test", las=2, t ="h")
  axis(1, at=(chr.lengths[,2] +  chr.lengths[,3]) / 2, paste0("Chr",rownames(chr.lengths)))
  abline(h=-log10(0.1 / N), col="orange") ; abline(h=-log10(0.01 / N), col="green")
  box()
  
  plot(x= map.autosomes[,"cumPos"],  y = -log10(counts[,"pMat"]), col = chr.cols, pch=19,xaxt='n', ylab="-log10(Pmat)", main = "MAT test", las=2, t ="h")
  axis(1, at=(chr.lengths[,2] +  chr.lengths[,3]) / 2, paste0("Chr",rownames(chr.lengths)))
  abline(h=-log10(0.1 / N), col="orange") ; abline(h=-log10(0.01 / N), col="green")
  box()
  
  plot(x= map.autosomes[,"cumPos"],  y = -log10(counts[,"pPoO"]), col = chr.cols, pch=19,xaxt='n', ylab="-log10(Ppofo)", main = "PofO test", las=2, t ="h")
  axis(1, at=(chr.lengths[,2] +  chr.lengths[,3]) / 2, paste0("Chr",rownames(chr.lengths)))
  abline(h=-log10(0.1 / N), col="orange") ; abline(h=-log10(0.01 / N), col="green")
  box()
  
  plot(x= map.autosomes[,"cumPos"],  y = -log10(counts[,"pPh"]), col = chr.cols, pch=19,xaxt='n', ylab="-log10(Phet)", main = "HET test", las=2, t ="h")
  axis(1, at=(chr.lengths[,2] +  chr.lengths[,3]) / 2, paste0("Chr",rownames(chr.lengths)))
  abline(h=-log10(0.1 / N), col="orange") ; abline(h=-log10(0.01 / N), col="green")
  if(!is.null(HWEdata)){
    plot(x= map.autosomes[,"cumPos"],  y = -log10(HWEdata[,"HWP"]), col = chr.cols, pch=19,xaxt='n', ylab="-log10(HW)", main = "HW equilibrium test", las=2, t ="h")
    axis(1, at=(chr.lengths[,2] +  chr.lengths[,3]) / 2, paste0("Chr",rownames(chr.lengths)))
    abline(h=-log10(0.1 / N), col="orange") ; abline(h=-log10(0.01 / N), col="green")
  }
  box()
}

HWEdata <- read.table("Analysis/HWEdata28.txt", sep = "\t")

HWEf2 <- HWEdata[,"HWPbh"] # ??
names(HWEf2) <- rownames(HWEdata)

create.plot(counts27)
create.plot(counts28, HWEdata)

map <- marker.annot[rownames(counts28), ]
map.autosomes <- map[which(map[,"Chr"] %in% 1:19),]

## Create the chromosome plot showing ATB and HWE

postscript("Figure_2.eps", horizontal = FALSE, paper = "special", width=12, height=8)

op <- par(mfrow=c(1,1))
op <- par(mar =c(4,4,2,1))

ymax <- max(as.numeric(map.autosomes[,"Pos"]))

plot(c(1,19), c(0, ymax+2000000), t = 'n', mgp=c(2.5,3.5,0), xlab="Chromosome", ylab="Position (Mbp)", xaxt='n', yaxt='n', main="Transmission bias from heterozygous parents", yaxs = "i")
axis(1, at=1:19, paste0(1:19), las=1, cex.axis=1)
axis(2, at=seq(0,ymax, 20000000), seq(0, ymax, 20000000) / 1000000, las=2, cex.axis=1)
for(x in 1:19){
  onChr <- rownames(map.autosomes[which(map.autosomes[,"Chr"] == as.character(x)),])
  colz <- as.numeric(-log10(counts28[onChr, "pPat"]) >  -log10(0.01/(70000*4))) + 1
  colz <- colz + as.numeric(as.character(counts28[onChr,"bfmi"]) == as.character(counts28[onChr,"patPref"]))
  colz[(-log10(counts28[onChr, "pPat"]) <=  -log10(0.01/(70000*4)))] <- 1
  
  colfunc <- c("white", "cornflowerblue", "orange")
  points(rep(x-0.15,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col=colfunc[colz], cex=1.8)

  colz <- as.numeric(-log10(counts28[onChr, "pMat"]) >  -log10(0.01/(70000*4))) + 1
  colz <- colz + as.numeric(as.character(counts28[onChr,"bfmi"]) == as.character(counts28[onChr,"matPref"]))
  colz[(-log10(counts28[onChr, "pMat"]) <=  -log10(0.01/(70000*4)))] <- 1

#  colfunc <- c("white", "gray30", "orange")
  points(rep(x+0.15,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col=colfunc[colz], cex=1.8)

  onC <- which(marker.annot[,"Chr"] == x)
  points(rep(x,length(onC)), as.numeric(marker.annot[onC,"Pos"]), pch="-",cex=1.5, col="gray")
  mcolz <- rep(1,length(onChr))
  mcolz[which(p.adjust(HWEf2[onChr]) < 0.01)] <- 2  
  points(rep(x,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col=c('black',"coral4")[mcolz],cex=1.5)
}
for(x in 1:19){
  mtext("♂", at = x - 0.22, side=1, cex=1.2)
  mtext("♀", at = x + 0.22, side=1, cex=1.2)
}
legend("topright", c("Suitable marker", "Unsuitable marker", "Marker out of HWE"), fill = c("black", "gray", "coral4"))
legend(x=16.8,y = 170000000, c("C57BL/6NCrl", "BFMI860-S12"), fill = c("cornflowerblue", "orange"))


dev.off()



## Create the chromosome plot for generation 28 and generation 27
map <- marker.annot[unique(rownames(counts27), rownames(counts28)), ]
map.autosomes <- map[which(map[,"Chr"] %in% 1:19),]


op <- par(mfrow=c(1,1))
ymax <- max(as.numeric(map.autosomes[,"Pos"]))

plot(c(1,19), c(0, ymax), t = 'n', xlab="Chromosome", ylab="Position (Mbp)", xaxt='n', yaxt='n', main="Transmission bias from heterozygous parents")
axis(1, at=1:19, paste0("Chr ", 1:19), las=2, cex.axis=0.7)
axis(2, at=seq(0,ymax, 20000000), seq(0, ymax, 20000000) / 1000000, las=2, cex.axis=0.9)
for(x in 1:19){
  onChr <- rownames(map.autosomes[which(map.autosomes[,"Chr"] == as.character(x)),])
  # generation 27
  points(rep(x - 0.15,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col="gray", cex=1.8)
  
  # generation 28
  colz <- as.numeric(-log10(counts28[onChr, "pPat"]) >  -log10(0.01/(70000*4))) + 1
  colz <- colz + as.numeric(-log10(counts28[onChr, "pMat"]) >  -log10(0.01/(70000*4)))
  points(rep(x + 0.15,length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col=colz, cex=1.8)
  
  points(rep(x+length(onChr)), as.numeric(map.autosomes[onChr,"Pos"]), pch="-", col="gray", cex=1.8)
}


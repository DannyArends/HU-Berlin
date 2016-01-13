# Analysis of Arabian horse SNP and performance data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2016
# first written Feb, 2015

toNumeric <- function(genotypes){
  numericG <- apply(genotypes, 1, function(x){
    geno <- table(unlist(strsplit(as.character(x),"")))
    a1 <- paste0(names(geno)[1],names(geno)[1])
    a2 <- paste0(sort(c(names(geno)[1],names(geno)[2])),collapse="")
    a3 <- paste0(names(geno)[2],names(geno)[2])
    ngeno <- rep(NA,length(x))
    ngeno[x == a1] <- 1
    ngeno[x == a2] <- 2
    ngeno[x == a3] <- 3
    return(ngeno)
  })
  rownames(numericG) <- colnames(genotypes)
  return(t(numericG))
}

toAB <- function(genotypes){
  numericG <- apply(genotypes, 1, function(x){
    geno <- table(unlist(strsplit(as.character(x),"")))
    a1 <- paste0(names(geno)[1],names(geno)[1])
    a2 <- paste0(sort(c(names(geno)[1],names(geno)[2])),collapse="")
    a3 <- paste0(names(geno)[2],names(geno)[2])
    ngeno <- rep(NA,length(x))
    ngeno[x == a1] <- "AA"
    ngeno[x == a2] <- "AB"
    ngeno[x == a3] <- "BB"
    return(ngeno)
  })
  rownames(numericG) <- colnames(genotypes)
  return(t(numericG))
}

setwd("E:/Horse/DNA/Kabadiner/")
# Read reference 'kabadiner' map and data
kabadinermap            <- read.table(file="input/cleaned_map.txt", sep = "\t")
kabadinergenotypes      <- read.table(file="input/cleaned_genotypes.txt", sep = "\t")
kabadinerphenotypes     <- read.table(file="input/cleaned_phenotypes.txt", sep = "\t")

# General horse chromosome info
chrInfo <- read.table("info/chrinfo.txt", sep="\t", header=TRUE)

setwd("E:/Horse/DNA/Equine60k/")
# Read data an pre-process
arabian <- read.table("input/arabianhorses.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=1)

phenotypes            <- arabian[1:28, ]                              # Row 1 till 29 contains phenotype data
genotypes             <- arabian[29:nrow(arabian), ]                  # Row 30 till the end contains genotype data
rownames(genotypes)   <- arabian[29:nrow(arabian), 1]
map                   <- genotypes[,c("Chromosome","Position")]       # Extract the map and split the phenotypes and genotypes
cat("Starting with", nrow(genotypes), "markers\n")
write.table(table(map[,"Chromosome"]), "output/AllMarkerTable.txt", sep="\t")
genotypes             <- genotypes[which(as.numeric(genotypes[,"GenTrain.Score"]) > 0.6),]    # Keep only high quality calls
map                   <- map[rownames(genotypes), ]

calledgeno  <- genotypes[,grep("GType", colnames(genotypes))]
genotypes   <- genotypes[,c(grep("Top.Alleles", colnames(genotypes)),grep("TOP", colnames(genotypes)))]
phenotypes  <- phenotypes[,grep("GType", colnames(phenotypes))]

colnames(phenotypes) <- gsub(".GType", "", colnames(phenotypes))
colnames(calledgeno) <- gsub(".GType", "", colnames(calledgeno))
colnames(genotypes)  <- gsub(".Top.Alleles", "", colnames(genotypes))
colnames(genotypes)  <- gsub(".TOP", "", colnames(genotypes))

### Data QC
cat("Left with", nrow(genotypes), "markers\n")
tables          <- apply(genotypes, 1, table)
nonInformative  <- which(unlist(lapply(tables,length)) < 2)             # Non informative, since we only have 1 genotype
genotypes       <- genotypes[-nonInformative, ]
calledgeno      <- calledgeno[-nonInformative, ]
map             <- map[-nonInformative, ]
cat("Left with", nrow(genotypes), "markers\n")

notDuplicated <- which(!duplicated(calledgeno))                         # Duplicated markers
genotypes     <- genotypes[notDuplicated, ]
map           <- map[notDuplicated, ]
cat("Left with", nrow(genotypes), "markers\n")

percMissing <- apply(apply(genotypes,1,is.na),2,sum) / ncol(genotypes)  # Call rates
keep          <- which(percMissing <= 0.05)
genotypes     <- genotypes[keep, ]
map           <- map[keep, ]
cat("Left with", nrow(genotypes), "markers\n")

allelefreq <- apply(genotypes, 1 , function(x){
  tbl <- table(unlist(lapply(x, strsplit, "")))
  min(tbl / sum(tbl))
})

keep <- which(allelefreq >= 0.05)
genotypes     <- genotypes[keep, ]
map           <- map[keep, ]
cat("Left with", nrow(genotypes), "markers\n")

write.table(map, file="input/cleaned_map.txt", sep = "\t")                                      # Save the clean map to disk
write.table(genotypes, file="input/cleaned_genotypes.txt", sep = "\t")                          # Save the clean genotypes to disk
write.table(toNumeric(genotypes), file="input/cleaned_numeric_genotypes.txt", sep = "\t")       # Encode the genotypes to be used for QTL mapping
write.table(phenotypes, file="input/cleaned_phenotypes.txt", sep = "\t")                        # Save the clean phenotypes to disk

# We now clean up the phenotypes manually after this step (fix stuff like "20 " -> "20") and only keep a single date for racing
phenotypes     <- read.table(file="input/cleaned_phenotypes_man.txt", sep = "\t", colClasses="character")

przewalski_ref   <- c("P0072", "P0078", "P0084")
write.table(kabadinergenotypes[, przewalski_ref], "input/cleaned_przewalski_genotypes.txt", sep="\t",quote=FALSE)

kabardian_ref    <- c("P3728", "P3421", "P3007", "P3542", "P3454", "P3509", "P3634", "P3420", "P3567", "P3418")
write.table(kabadinergenotypes[, kabardian_ref], "input/cleaned_kabardian_genotypes.txt", sep="\t",quote=FALSE)

thoroughbred_ref <- c("P3708", "P0031", "P1294", "P3706", "P3751", "P3752")
write.table(kabadinergenotypes[, thoroughbred_ref], "input/cleaned_thoroughbred_genotypes.txt", sep="\t",quote=FALSE)

# Which markers match
kabadinermap <- kabadinermap[which(rownames(kabadinermap) %in% rownames(map)),]

# Sort the map and the genotypes of the arabian horses so that they match the kabadinermap
map                 <- map[rownames(kabadinermap),]
genotypes           <- genotypes[rownames(kabadinermap), ]
kabadinergenotypes  <- kabadinergenotypes[rownames(kabadinermap), ]
cat("Left with", nrow(genotypes), "markers\n")
write.table(table(map[,"Chromosome"]), "output/GoodMarkerTable.txt", sep="\t")

genotypesnref <- cbind(genotypes, kabadinergenotypes[,przewalski_ref], kabadinergenotypes[,kabardian_ref], kabadinergenotypes[,thoroughbred_ref])
cat("Shared:", nrow(genotypesnref), "markers\n")

for(x in unique(map[,"Chromosome"])){
  cat(x, length(which(map[,"Chromosome"] == x)),"\n")
}

numGeno <- t(toNumeric(genotypesnref))

if(!file.exists("input/cleaned_genotypes_structure.txt")){
  # Write out the data for STRUCTURE
  structGeno <- NULL #matrix(NA, nrow(numGeno) * 2, ncol(numGeno))
  for(x in 1:nrow(numGeno)){
    gg <- rbind(rep(NA, ncol(numGeno)), rep(NA, ncol(numGeno)))
    a1 <- which(numGeno[x,] == 1)
    a2 <- which(numGeno[x,] == 2)
    a3 <- which(numGeno[x,] == 3)
    gg[1, a1] <- 0; gg[2, a1] <- 0  # Pretty inefficient, but it will do the trick
    gg[1, a2] <- 0; gg[2, a2] <- 1
    gg[1, a3] <- 1; gg[2, a3] <- 1
    gg[is.na(gg)] <- 9
    structGeno <- rbind(structGeno, gg)
  }

  rownames(structGeno) <- unlist(lapply(rownames(numGeno), rep, 2))
  colnames(structGeno) <- colnames(numGeno)
  write.table(structGeno, file="input/cleaned_genotypes_structure.txt", sep = "\t")    # Save the clean genotypes to disk
}

## Some basic plots of all the individuals relatedness
dendrogram <- as.dendrogram(hclust(dist(t(toNumeric(genotypesnref)), method = "manhattan")))         #TODO: perhaps add the reference horse

strains <- c(as.character(phenotypes["Strain", ]), rep("P", length(przewalski_ref)), rep("Kab", length(kabardian_ref)), rep("Tho", length(thoroughbred_ref)))
names(strains) <- c(colnames(phenotypes), przewalski_ref, kabardian_ref, thoroughbred_ref)

# Create colors
cols <- c("red", "blue", "orange", "black", "purple", "brown")
names(cols) <- c("K", "S", "H", "P", "Kab", "Tho")

labelCol <- function(x) {
  if (is.leaf(x)) {
    hclass <- strains[attr(x, "label")]             # Fetch the class label
    hcol <- cols[hclass]                            # Determine color of the label
    cat(attr(x, "label"), hclass, hcol, "\n")
    attr(x, "nodePar") <- list(lab.col=hcol)
  }
  return(x)
}
# Add the colors to the dendrogram
dendrogram.col <- dendrapply(dendrogram, labelCol)
plot(dendrogram.col, main = "")

# Wright statistics of every marker (overall and per group)
if(!file.exists("output/WrightFstatistics.txt")) {
  Fvalues <- rep(NA, nrow(genotypes))             # Use only the genotypes of the Arabian horses
  names(Fvalues) <- rownames(genotypes)

  # Wright statistics
  doFtest <- function(genotypes, group) {
    for(m in rownames(genotypes)){
      N <- sum(!is.na(genotypes[m, group]))                                       # How many non NA individuals
      splitgeno <- lapply(genotypes[m, group], function(x) {
        strsplit(as.character(x), "")}                                            # Split "AA" -> "A", "A"
      )
      genotype <- table(unlist(splitgeno))                                        # Number of A / T or C / G
      genotype <- genotype / sum(genotype)                                        # Genotype % of A / T
      expected <- N * (2*genotype[1]*genotype[2])                                 # Expected heterozygous in HW
      isH <- unlist(lapply(splitgeno, function(y){y[[1]][1] != y[[1]][2]}))
      F = 1 - (sum(as.numeric(isH), na.rm=TRUE) / expected)                       # F-test
      cat(m, expected, sum(as.numeric(isH), na.rm=TRUE), F, "\n")
      Fvalues[m] <- F
    }
    return(Fvalues)
  }

  Fstatistics <- cbind(
    populationF = doFtest(genotypes, colnames(phenotypes)),
    HamdaniF = doFtest(genotypes, colnames(phenotypes[, phenotypes["Strain",] == "H"])),
    SaklawiF = doFtest(genotypes, colnames(phenotypes[, phenotypes["Strain",] == "S"])),
    KahlawiF = doFtest(genotypes, colnames(phenotypes[, phenotypes["Strain",] == "K"]))
  )
  write.table(Fstatistics, "output/WrightFstatistics.txt", sep="\t")
}else{
  Fstatistics <- read.table("output/WrightFstatistics.txt", sep="\t")
}

### Chromosome plot to show the location of SNPs
chromosomes  <- as.character(c(1:31, "X", "Y", "MT"))

plot(y=c(0, max(chrInfo[,2])), x=c(1,nrow(chrInfo)), t='n', main="", yaxt="n", xlab="Chromosome", ylab="Length (Mb)", xaxt="n")
abline(h=seq(0, max(chrInfo[,2]), 10000000), col = "lightgray", lty = "dotted")

cnt <- 1
aa <- apply(chrInfo,1,function(x) {
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1,lwd=3)
  cnt <<- cnt + 1
})

aa <- apply(cbind(map, F = Fstatistics), 1 ,function(x) {
  xloc <- match(as.character(x["Chromosome"]), chromosomes); yloc <- as.numeric(x["Position"])
  col <- "gray"
  if(as.numeric(x["F.populationF"]) < -0.25) col <- "lightblue"
  if(as.numeric(x["F.populationF"]) < -0.50) col <- "blue"
  if(as.numeric(x["F.populationF"]) >  0.25) col <- "yellow"
  if(as.numeric(x["F.populationF"]) >  0.50) col <- "red"
  points(x=xloc, y=yloc, pch="-",cex=1.6, col=col)
})

axis(1, chromosomes, at=c(1:nrow(chrInfo)), las=1)
axis(2, seq(0, max(chrInfo[,2]), 10000000)/1000000, at=seq(0, max(chrInfo[,2]), 10000000), cex.axis=0.7)

# Pairwise population D and Fst values using StAMPP
library(StAMPP)
abGeno <- t(toAB(genotypesnref))
stammpinput <- abGeno
#strains <- as.character(unlist(c(phenotypes["Strain",rownames(stammpinput)[1:48]], "P", "P", "P")))
stammpinput <- data.frame(cbind(rownames(stammpinput), strains, 2, "BiA", stammpinput))
colnames(stammpinput)[1:4] <- c("Sample", "Pop", "Ploidy", "Format")

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values
#    S         K        H        P
#S 0.000000 0.018056 0.022127 0.162171
#K 0.018056 0.000000 0.023948 0.163087
#H 0.022127 0.023948 0.000000 0.172120
#P 0.162171 0.163087 0.172120 0.000000

stammpinput.fst <- stamppFst(stammpinput.freq, 1000, 95, 4) # Population Fst values
#            S           K         H  P
#S          NA          NA        NA NA
#K 0.002999697          NA        NA NA
#H 0.008421684 0.009359349        NA NA
#P 0.196936652 0.196610481 0.2122905 NA

stammp.D.ind <- stamppNeisD(stammpinput.freq, FALSE)    # Distance between individuals
stamppAmova(stammp.D.ind, stammpinput.freq, 100)        # Calculate AMOVA

## Structure results (run Structure)
structuredir <- "E:/Horse/DNA/Equine60k/STRUCTURE"
projectname <- "ArabianHorses"
paramsetname <- "5000-1000"

loc <- paste0(structuredir, "/", projectname, "/", paramsetname, "/Results")
results <- dir(loc)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)       # Returns string w/o leading or trailing whitespace

analyzeStructure <- function(stmatrix, confidence = 0.0){
  stmatrix <- stmatrix[which(apply(stmatrix,1,max) > confidence), ]
  mp <- apply(stmatrix, 1, which.max)
  mp <- cbind(mp, Strain = NA)

  ii <- rownames(mp)[which(rownames(mp) %in% colnames(phenotypes))]
  mp[ii, "Strain"] <- as.character(phenotypes["Strain", ii])
  mp[is.na(mp)] <- "P"

  counts <- matrix(0, ncol(stmatrix), 4, dimnames=list(1:ncol(stmatrix), c("H","S","K", "P")))
  for(x in 1:nrow(mp)){
    counts[as.numeric(mp[x,1]), mp[x,2]] <- counts[mp[x,1], mp[x,2]] + 1
  }
  print(counts)
  return(counts)
}

plotStructure <- function(stmatrix, doSort = FALSE){
  if(doSort){
    ordering <- NULL
    breaks <- 0
    stmatrix.copy <- stmatrix
    for(x in 1:ncol(stmatrix.copy)){
      sorted <- sort(stmatrix.copy[,x], decreasing = TRUE)
      samples <- names(sorted[sorted >= 0.5])
      ordering <- c(ordering, samples)
      stmatrix.copy <- stmatrix.copy[-which(rownames(stmatrix.copy) %in% samples),]
      breaks <- c(breaks, breaks[length(breaks)] + length(samples))
      cat(breaks, samples,"\n")
    }
    #breaks <- c(breaks, nrow(stmatrix))
    ordering <- c(ordering, rownames(stmatrix.copy))
    stmatrix <- stmatrix[ordering,]
  }
  
  plot(c(1,nrow(stmatrix)), c(-0.15, 1), t = 'n', xaxt='n', xlab = "Individual", ylab = "Cluster membership (%)", yaxt='n', main=paste0("STRUCTURE, clusters = ", ncol(stmatrix)))
  dsum <- rep(0, nrow(stmatrix))
  mcol <- 2
  apply(stmatrix, 2, function(x){
    for(i in 1:length(x)){
      rect(i-0.4, dsum[i], i+0.4, dsum[i] + x[i], col=mcol)
    }
    dsum <<- dsum + x
    mcol <<- mcol + 1
    #cat(dsum, "\n")
  })
  sahriastrain <- rep("P", length(rownames(stmatrix)))
  names(sahriastrain) <- rownames(stmatrix)
  sahriastrain[names(phenotypes["Strain",])] <- as.character(phenotypes["Strain",])
  text(x = 1:nrow(stmatrix), y = rep(-0.05, nrow(stmatrix)), sahriastrain)
  mids <- diff(breaks) / 2 + 0.5
  for(x in 1:length(mids)) mids[x] <- mids[x] + breaks[x]
  text(x = mids, y = rep(-0.1, length(mids)-1), paste0("Group ", 1:length(mids-1)), cex=0.8) 
  axis(1, at=1:nrow(stmatrix), rownames(stmatrix), las = 2, cex.axis = 0.8)
  axis(2, at=seq(0, 1, 0.1), seq(0, 100, 10), las = 2, cex.axis = 0.8)
  abline(v = breaks + 0.5,lwd=0.5, lty=2)
  return(breaks)
}

for(analysis in paste0(loc, "/", results)){
  structuredata <- readLines(analysis, warn = FALSE)
  bt <- which(grepl("Estimated Ln Prob of Data", structuredata))                        # Start of info
  st <- which(grepl("Inferred ancestry of individuals:", structuredata))+2              # Start of data table
  et <- which(grepl("Estimated Allele Frequencies in each cluster", structuredata))-3   # End of data table
  
  datalines <- strsplit(trim(unlist(lapply(strsplit(structuredata[st:et], ":"), "[", 2))), " ")
  stmatrix <- NULL
  for(l in datalines){
    stmatrix <- rbind(stmatrix, as.numeric(strsplit(l, " ")))
  }
  rownames(stmatrix) <- unlist(lapply(strsplit(structuredata[st:et], "\""),"[", 2))
  cat("--", analysis, "--\n")
  cat(structuredata[bt:(st-5)],sep="\n")        # Print some structure information
  plotStructure(stmatrix, TRUE)
  cat("--All individuals--\n")
  counts <- analyzeStructure(stmatrix)          # Compare how good the structure model fits with the breeders perspective
  cat("--High purity (> 0.8)--\n")
  counts <- analyzeStructure(stmatrix, 0.8)     # Compare how good the pure structure model fits with the breeders perspective
  cat("----\n")
  if(ncol(stmatrix) == 4) break                 # We did more but I want too see if Saria's hypothesis is correct
}

# Analyse the different covariates for the different phenotypes
sex     <- as.factor(unlist(phenotypes["Sex",]))
strain  <- as.factor(unlist(phenotypes["Strain",]))
ageAtM  <- as.numeric(unlist(phenotypes["D.Measure",])) - as.numeric(unlist(phenotypes["D.Birth",]))
ageAtR  <- as.numeric(unlist(phenotypes["D.Racing",])) - as.numeric(unlist(phenotypes["D.Birth",]))

classicalpheno <- c("WH", "CW", "CH", "NG", "TG", "ChG", "ChD", "ChW", "BLL", "BL", "FCL", "HCL")
racepheno <- c("Distance (km)", "Speed(km\\hr)")

pClassical <- matrix(NA, length(classicalpheno), 3, dimnames = list(classicalpheno, c("Sex", "Age", "Strain")))
for(ph in classicalpheno) {
  model <- anova(lm(as.numeric(phenotypes[ph,]) ~ sex + ageAtM + strain))
  pClassical[ph,] <- model[[5]][1:3]
}

pRace <- matrix(NA, length(racepheno), 3, dimnames=list(racepheno, c("Sex", "Age", "Strain")))
for(ph in racepheno) {
  model <- anova(lm(as.numeric(phenotypes[ph,]) ~ sex + ageAtR + strain))
  pRace[ph,] <- model[[5]][1:3]
}

write.table(rbind(pClassical, pRace), "output/covariates.txt", sep="\t")

# Haplotype analysis

pred <- as.factor(as.character(phenotypes["mtDNA HT",]))

pClassical <- matrix(NA, length(classicalpheno), 3, dimnames = list(classicalpheno, c("Sex", "Age", "mtDNA HT")))
for(ph in classicalpheno) {
  model <- anova(lm(as.numeric(phenotypes[ph,]) ~ sex + ageAtM + pred))
  pClassical[ph,] <- model[[5]][1:3]
}

pRace <- matrix(NA, length(racepheno), 3, dimnames=list(racepheno, c("Sex", "Age", "mtDNA HT")))
for(ph in racepheno) {
  model <- anova(lm(as.numeric(phenotypes[ph,]) ~ sex + ageAtR + pred))
  pRace[ph,] <- model[[5]][1:3]
}

write.table(rbind(pClassical, pRace), "output/mtDNA_haplotypes.txt", sep="\t")

# GWAS / QTL analysis

# Calculate the phenotypes after correction
phenoC <- matrix(NA, length(c(classicalpheno, racepheno)), ncol(phenotypes), dimnames = list(c(classicalpheno, racepheno), colnames(phenotypes)))
for(ph in classicalpheno) {
  model <- lm(as.numeric(phenotypes[ph,]) ~ sex + ageAtM + strain + pred)
  pCorrected <- model$residuals + model$coefficients["(Intercept)"]
  phenoC[ph, as.numeric(names(pCorrected))] <- pCorrected
}
for(ph in racepheno) {
  model <- lm(as.numeric(phenotypes[ph,]) ~ sex + ageAtR + strain + pred)
  pCorrected <- model$residuals + model$coefficients["(Intercept)"]
  phenoC[ph, as.numeric(names(pCorrected))] <- pCorrected
}

# Do the GWAS using a single QTL model *Do not correct for race of horse
pvalues <- matrix(NA, length(c(classicalpheno, racepheno)), nrow(genotypes), dimnames = list(c(classicalpheno, racepheno), rownames(genotypes)))
for(phe in rownames(phenoC)) {
  cat("Computing GWAS results for:", phe, "\n")
  pvalues[phe, ] <- apply(genotypes, 1, function(marker) {
    tryCatch(res <- anova(lm(as.numeric(phenoC[phe,]) ~ as.factor(marker)))[[5]][1], error = function(e){ res <<- NA })
    return(res)
  })
}
pvalues <- t(pvalues)
write.table(pvalues, "output/pvaluesGWAS.txt", sep="\t")

setwd("E:/Horse/DNA/Equine60k/")
pvalues <- read.csv("output/pvaluesGWAS.txt", sep="\t")

chrcols <- 1+ (as.numeric(as.factor(map[,"Chromosome"])) %% 2)

for(phe in colnames(pvalues)) {
  ii <- which(pvalues[, phe] < (0.05/nrow(pvalues)))
  if(length(ii) > 0){
    cat(phe, rownames(pvalues)[ii],"\n")
  }
}

for(phe in colnames(pvalues)) {
  png(paste0("output/GWAS_", gsub("\\\\hr", "pH", phe), ".png"), width=1024, height=600)
    plot(x = c(1, nrow(pvalues)), y = c(0,10), t='n', main=phe)
    points(-log10(pvalues[,phe]), pch = 19, cex = 0.5, col=chrcols)
    abline(h = -log10(0.1/nrow(pvalues)), col="orange")
    abline(h = -log10(0.05/nrow(pvalues)), col="gold")
    abline(h = -log10(0.01/nrow(pvalues)), col="green")
  dev.off()
}

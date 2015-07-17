# Additional SNPs
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2015
# first written Jun, 2015

source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")

kasp   <- read.table("Analysis/kasp+licor.txt", sep="\t", check.names=FALSE, colClasses="character", skip=2, header=TRUE, na.strings=c("","-", "?","NA", "CG?", "CT?", "CT?", "AG?", "KL?"))

map <- t(kasp[c(1,3), -c(1:6)])
colnames(map) <- c("Chr", "Mb_NCBI38")

kasp <- kasp[-c(1:3), -c(2:6)]
rownames(kasp) <- kasp[,1]
kasp <- kasp[,-1]

goodM <- colnames(kasp)[ which(apply(kasp,2,function(x){return(length(which(is.na(x)))) }) < 400) ] # colnames of good markers

genotypes <- t(kasp[,goodM])   # Individuals in the columns, SNPs in the rows
map <- map[goodM,]        # SNPs in the rows, Chromosome and Position (NCBI38)

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("33310233", "6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                         # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]; length(F2)                                     # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]; length(F1)                                     # The F1 individuals
P <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]; length(P)                                       # The P individuals

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
birthmonth <- unlist(lapply(strsplit(as.character(phenotypes[,"W.dat"]),".", fixed=TRUE),"[",2))
phenotypes <- cbind(phenotypes, Birthmonth = birthmonth)                                                      # Add the birth month column to the matrix

phenotypes <- cbind(phenotypes, mri42d_fatDlean = phenotypes[,"mri42d_fat"] / phenotypes[,"mri42d_lean"])     # Fat / Lean day 42
phenotypes <- cbind(phenotypes, mri56d_fatDlean = phenotypes[,"mri56d_fat"] / phenotypes[,"mri56d_lean"])     # Fat / Lean day 56
phenotypes <- cbind(phenotypes, mri70d_fatDlean = phenotypes[,"mri70d_fat"] / phenotypes[,"mri70d_lean"])     # Fat / Lean day 70

genotypesN <- apply(genotypes, 1, function(x){
  cat(x["BFMI"], x["B6N"],"\n")
  a <- rep("",length(x))
  a[!is.na(x)] <- "H"
  a[as.character(x) == as.character(x["BFMI"])] <- "BFMI"
  a[as.character(x) == as.character(x["B6N"])] <- "B6N"
  return(a)
})
rownames(genotypesN) <- colnames(genotypes)

genotypes[,c("6661385","6661114","6662242","6661117","BFMI","B6N")]
t(genotypesN)[,c("6661385","6661114","6662242","6661117","BFMI","B6N")]


phenotypes <- phenotypes[F2,]
genotypes <- t(genotypesN)[,F2]

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
for(phe in phenonames){
  pvalues <- NULL
  for(x in 1:nrow(genotypes)){
    ind            <- names(genotypes[x,!is.na(genotypes[x,])])                                                            # Which individuals have genotype data

    subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure (factor)
    littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter  (linear effect)
    litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter    (factor)
    season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born    (factor)
    genotype       <- as.factor(t(genotypes[x,!is.na(genotypes[x,])]))                                                        # The genotype under investigation  (factor)

    phenotype      <- phenotypes[ind, phe]                      #/ phenotypes[ind, paste0("mri",pheno.col,"_lean")]     # Response: Fat / Lean
    tryCatch(res <- anova(lm(phenotype ~ subfamily + littersize + litternumber + season + genotype))[[5]], error = function(e){ res <<- rep(NA, 5) })
    cat(rownames(genotypes)[x], map[rownames(genotypes)[x], "Mb_NCBI38"], round(-log10(res),2),"\n")
    pvalues <- rbind(pvalues, -log10(res[-length(res)]))
  }
  rownames(pvalues) <- rownames(genotypes)
  colnames(pvalues) <- c("subfamily", "l_size", "l_number", "season", "marker")
  write.table(pvalues, paste0("Analysis/Kasp/QTL_",phe,"_kasp.txt"),sep="\t")
}

setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga
mapO <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypesO   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes
map <- rbind(map, mapO[,1:2])

genotypesOG <- apply(genotypesO, 1, function(x){
  cat(x["BFMI860-12 (V2)"], x["B6N"],"\n")
  a <- rep("",length(x))
  a[as.character(x) == as.character(x["BFMI860-12 (V2)"])] <- "BFMI"
  a[as.character(x) != as.character(x["BFMI860-12 (V2)"])] <- "B6N"
  a[as.character(x) == "H"] <- "H"
  return(a)
})

rownames(genotypesOG) <- colnames(genotypesO)
genotypesOG <- t(genotypesOG)

genotypes[,c("6661385","6661114","6662242","6661117")]


setwd("E:/Mouse/ClassicalPhenotypes/AIL")
qtl42FatO  <- read.table("Analysis/qtls_mri70d_fat.txt",   sep="\t")
qtl42Fat   <- read.table("Analysis/kasp/QTL_mri70d_fat_kasp.txt",   sep="\t")
qtl42Fat   <- rbind(qtl42Fat, qtl42FatO)


submap <- map
submap <- map[rownames(qtl42Fat[qtl42Fat[,"marker"] > 1,]), ]
submap <- submap[submap[,"Chr"] == 3,]                                    #  & 
submap <- submap[rownames(submap) %in% rownames(qtl42Fat),]


submap <- submap[as.numeric(as.character(submap[,"Mb_NCBI38"])) > 32500000 &  as.numeric(as.character(submap[,"Mb_NCBI38"])) < 40000000,]
submap <- submap[order(as.numeric(as.character(submap[,2]))),]

colz <- rep("black", nrow(submap))
colz[grep("KM", rownames(submap))] <- "red"

#op <- par(mfrow=c(1,3))
png("Analysis/Region_Chr3_Before.png")
op <- par(cex=1.5)
plot(c(32500000, 40000000), y = c(0, 65), t = 'n', ylab = "LOD", xlab = "Chromosome 3: 30 Mb - 45 Mb", xaxt='n', las = 2, main= "QTL  + Lab markers")
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = qtl42Fat[rownames(submap), "marker"], t = 'l', col = "red",lwd=3)
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = rep(-1.3, nrow(submap)), pch="|", col = colz, cex=0.5)
axis(1, at=seq(32500000, 40000000, 2500000), seq(32500000, 40000000, 2500000) / 1000000)
dev.off()

submapO <- submap[-grep("KM", rownames(submap)),]
png("Analysis/Region_Chr3_After.png")
op <- par(cex=1.5)
plot(c(32500000, 40000000), y = c(0, 65), t = 'n', ylab = "LOD", xlab = "Chromosome 3: 30 Mb - 45 Mb", xaxt='n', las = 2, main= "MegaMuga")
points(x = as.numeric(as.character(submapO[,"Mb_NCBI38"])), y = qtl42Fat[rownames(submapO), "marker"], t = 'l', col = "red",lwd=3)
points(x = as.numeric(as.character(submapO[,"Mb_NCBI38"])), y = rep(-1.3, nrow(submapO)), pch="|", col = "black", cex=0.5)
axis(1, at=seq(32500000, 40000000, 2500000), seq(32500000, 40000000, 2500000) / 1000000)
dev.off()

genotypesOG <- genotypesOG[,which(colnames(genotypesOG) %in% colnames(genotypes))]

genotypesOG <- genotypesOG[,match(colnames(genotypes), colnames(genotypesOG))]

#colnames(genotypesOG[,match(colnames(genotypes), colnames(genotypesOG))])[1:10]
#colnames(genotypes)[1:10]



combinedGeno <- rbind(genotypesOG,genotypes)

genotypes[,c("6661385","6661114","6662242","6661117")]
combinedGeno["KM-023",c("6661385","6661114","6662242","6661117")]

cat("Missing before:", length(which(combinedGeno == "")),"\n")

combinedGeno["UNC5009661", "6661069"]

combinedGeno <- combinedGeno[rownames(submap),]



x <- which(rownames(combinedGeno) == "UNC5009661")
y <- "6661069"

for(y in 1:ncol(combinedGeno)){
  for(x in 2:(nrow(combinedGeno)-1)){
    if(combinedGeno[x,y] == ""){
      #cat(x,y,":", combinedGeno[(x-1), y], combinedGeno[x,y],combinedGeno[(x+1),y],"\n")
      if(combinedGeno[(x-1), y] == combinedGeno[(x+1),y]){
        combinedGeno[x,y] <- combinedGeno[(x-1), y]
      }
    }
  }
}
cat("Missing after:", length(which(combinedGeno == "")),"\n")


setwd("E:/Mouse/ClassicalPhenotypes/AIL")
for(phe in phenonames){
  pvalues <- NULL
  for(x in 1:nrow(combinedGeno)){
    ind            <- names(combinedGeno[x,!is.na(combinedGeno[x,])])                                                            # Which individuals have genotype data

    subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure (factor)
    littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter  (linear effect)
    litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter    (factor)
    season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born    (factor)
    genotype       <- as.factor(t(combinedGeno[x,!is.na(combinedGeno[x,])]))                                                        # The genotype under investigation  (factor)

    phenotype      <- phenotypes[ind, phe]                      #/ phenotypes[ind, paste0("mri",pheno.col,"_lean")]     # Response: Fat / Lean
    tryCatch(res <- anova(lm(phenotype ~ subfamily + littersize + litternumber + season + genotype))[[5]], error = function(e){ res <<- rep(NA, 5) })
    cat(rownames(combinedGeno)[x], map[rownames(combinedGeno)[x], "Mb_NCBI38"], round(-log10(res),2),"\n")
    pvalues <- rbind(pvalues, -log10(res[-length(res)]))
  }
  rownames(pvalues) <- rownames(combinedGeno)
  colnames(pvalues) <- c("subfamily", "l_size", "l_number", "season", "marker")
  write.table(pvalues, paste0("Analysis/Kasp/QTLs_",phe,"_region.txt"),sep="\t")
}

qtl42Fat   <- read.table("Analysis/kasp/QTLs_mri42d_fat_region.txt",   sep="\t")
qtl70Fat   <- read.table("Analysis/kasp/QTLs_mri70d_fat_region.txt",   sep="\t")
qtl42   <- read.table("Analysis/kasp/QTLs_d42_region.txt",   sep="\t")
qtl70   <- read.table("Analysis/kasp/QTLs_d70_region.txt",   sep="\t")

#op <- par(mfrow=c(1,2))
png("Analysis/Region_Chr3_AfterSmoothed.png")
op <- par(cex=1.5)
plot(c(32500000, 40000000), y = c(0, 65), t = 'n', ylab = "LOD", xlab = "Chromosome 3: 30 Mb - 45 Mb", xaxt='n', las = 2, main= "MegaMuga + Lab markers + Smoothing")
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = qtl42Fat[rownames(submap), "marker"], t = 'l', col = "red",lwd=3)
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = rep(-1.3, nrow(submap)), pch="|", col = colz, cex=0.5)
axis(1, at=seq(32500000, 40000000, 2500000), seq(32500000, 40000000, 2500000) / 1000000)
dev.off()

# Debug print some data to sort out if everything is still equal to the RAW data
#
#combinedGeno["KM-023",c("6661385","6661114","6662242","6661117")]
#combinedGeno[c("UNC5009661","KM-023", "KM-025"),c("6661068")]
#genotypes[c("KM-023", "KM-025"),c("6661068")]
#combinedGeno["UNC5009661", "6661069"]

combinedGP <- rbind(phenotypes[colnames(combinedGeno),"Vater"], 
                    phenotypes[colnames(combinedGeno),"d42"], phenotypes[colnames(combinedGeno),"d70"] , 
                    phenotypes[colnames(combinedGeno),"mri42d_fat"], phenotypes[colnames(combinedGeno),"mri70d_fat"] , 
                    phenotypes[colnames(combinedGeno),"mri42d_lean"], phenotypes[colnames(combinedGeno),"mri70d_lean"] , 
                    combinedGeno[rownames(submap),])
#combinedGP <- rbind(phenotypes[colnames(combinedGeno),"d42"] , combinedGeno[rownames(submap),])

combinedGP["UNC5009661", "6661069"]
combinedGPL <- cbind( map[rownames(combinedGP),], qtl42Fat[rownames(combinedGP), "marker"], qtl70Fat[rownames(combinedGP), "marker"], qtl42[rownames(combinedGP), "marker"], qtl70[rownames(combinedGP), "marker"], combinedGP)


combinedGPL["UNC5009661", "6661069"]

rownames(combinedGPL)[1:7] <- c("Vater", "d42","d70", "mri42d_fat", "mri70d_fat", "mri42d_lean", "mri70d_lean")
colnames(combinedGPL)[3:6] <- c("LOD_mri_42d_fat", "LOD_mri_70d_fat", "LOD_42d", "LOD_70d")
write.table(combinedGPL, "Analysis/Chromosome3region_MarkersLodPheno.txt",sep="\t")
combinedGPL <- read.table("Analysis/Chromosome3region_MarkersLodPheno.txt", sep = "\t", colClasses="character",check.names=FALSE)

setwd("E:/Mouse/DNA/MegaMuga/")

LAB <- read.table("ManualMarkersLab.txt",sep='\t', header = TRUE, na.strings=c("","NA","Undetermined", " gg", "gg", " cc", "cc", "tt","aa"), colClasses="character")
LAB <- LAB[,c("TierID","Ccna2","Bbs7","M313")]
LAB[,1] <- paste0("666", LAB[,1])
LAB <- LAB[which(LAB[,1] %in% F2),]
rownames(LAB) <- LAB[,1]
LAB <- LAB[,-1]

locations <- c(36567121, 36575562, 36824165)

LAB[which(LAB[,1] == "GG"),1] <- "BFMI";LAB[which(LAB[,1] == "AG"),1] <- "H";LAB[which(LAB[,1] == "AA"),1] <- "B6N"
LAB[which(LAB[,2] == "GG"),2] <- "BFMI";LAB[which(LAB[,2] == "AG"),2] <- "H";LAB[which(LAB[,2] == "AA"),2] <- "B6N"
LAB[which(LAB[,3] == "CC"),3] <- "BFMI";LAB[which(LAB[,3] == "CT"),3] <- "H";LAB[which(LAB[,3] == "TT"),3] <- "B6N"
LAB <- t(LAB)
LAB <- LAB[,match(colnames(combinedGP), colnames(LAB))]

LAB <- cbind(3, locations, NA, NA, NA, NA, LAB)
colnames(LAB)[1:6] <- c("Chr", "Mb_NCBI38", "LOD_mri_42d_fat", "LOD_mri_70d_fat", "LOD_42d", "LOD_70d")

combinedALL <- rbind(LAB,combinedGPL)

PheSection <- combinedALL[4:10,]
GenSection <- combinedALL[-c(4:10),]

idx <- sort(as.numeric(as.character(GenSection[,2])), index.return=TRUE)$ix
GenSection <- GenSection[idx,]

combinedALLordered <- rbind(PheSection,GenSection)
combinedALLordered[,3] <- round(as.numeric(as.character(combinedALLordered[,3])),2)
combinedALLordered[,4] <- round(as.numeric(as.character(combinedALLordered[,4])),2)
combinedALLordered[,5] <- round(as.numeric(as.character(combinedALLordered[,5])),2)
combinedALLordered[,6] <- round(as.numeric(as.character(combinedALLordered[,6])),2)
setwd("E:/Mouse/ClassicalPhenotypes/AIL")
write.table(combinedALLordered, "Analysis/Chromosome3region_MarkersLodPheno.txt",sep="\t")

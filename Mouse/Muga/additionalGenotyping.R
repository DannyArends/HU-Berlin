# Additional SNPs
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2015
# first written Jun, 2015

source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")

kasp   <- read.table("Analysis/kasp.txt", sep="\t", check.names=FALSE, colClasses="character", skip=2, header=TRUE, na.strings=c("","-","NA", "CG?", "CT?", "CT?", "AG?"))

map <- t(kasp[c(1,3), -c(1:6)])
colnames(map) <- c("Chr", "Mb_NCBI38")

kasp <- kasp[-c(1:3), -c(2:6)]
rownames(kasp) <- kasp[,1]
kasp <- kasp[,-1]

goodM <- colnames(kasp)[ which(apply(kasp,2,function(x){return(length(which(is.na(x)))) }) < 25) ] # colnames of good markers

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

phenotypes <- phenotypes[F2,]
genotypes <- genotypes[,F2]

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
map <- rbind(map, mapO[,1:2])

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
qtl42FatO  <- read.table("Analysis/qtls_mri42d_fat.txt",   sep="\t")
qtl42Fat   <- read.table("Analysis/kasp/QTL_mri42d_fat_kasp.txt",   sep="\t")
qtl42Fat   <- rbind(qtl42Fat, qtl42FatO)


submap <- map
submap <- map[rownames(qtl42Fat[qtl42Fat[,"marker"] > 1,]), ]
submap <- submap[submap[,"Chr"] == 3,]                                    #  & 
submap <- submap[rownames(submap) %in% rownames(qtl42Fat),]


submap <- submap[as.numeric(as.character(submap[,"Mb_NCBI38"])) > 32500000 &  as.numeric(as.character(submap[,"Mb_NCBI38"])) < 40000000,]
submap <- submap[order(as.numeric(as.character(submap[,2]))),]

colz <- rep("black", nrow(submap))
colz[grep("KM", rownames(submap))] <- "red"

op <- par(mfrow=c(1,2))
plot(c(32500000, 40000000), y = c(0, 65), t = 'n', ylab = "LOD", xlab = "Chromosome 3: 30 Mb - 45 Mb", xaxt='n', las = 2, main= "After")
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = qtl42Fat[rownames(submap), "marker"], t = 'l', col = "red")
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = rep(-1.3, nrow(submap)), pch="|", col = colz, cex=0.5)
axis(1, at=seq(32500000, 40000000, 2500000), seq(32500000, 40000000, 2500000) / 1000000)

submap <- submap[-grep("KM", rownames(submap)),]

plot(c(32500000, 40000000), y = c(0, 65), t = 'n', ylab = "LOD", xlab = "Chromosome 3: 30 Mb - 45 Mb", xaxt='n', las = 2, main= "Before")
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = qtl42Fat[rownames(submap), "marker"], t = 'l', col = "red")
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = rep(-1.3, nrow(submap)), pch="|", col = "black", cex=0.5)
axis(1, at=seq(32500000, 40000000, 2500000), seq(32500000, 40000000, 2500000) / 1000000)








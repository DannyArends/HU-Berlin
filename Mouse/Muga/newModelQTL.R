#
# New QTL mapping of BFMI, using a slightly different population structure correction
#

source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")                                                                                           # Read in the data from the Mega Muga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE, colClasses="character")                  # Normal A, H, B genotypes

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes  <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)                      
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                                 # Add the season column to the matrix
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                                     # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                                        # Remove the faulty individuals from the genotypes
                                  
F2  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                            # The F2 individuals
F1  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]                                                            # The F1 individuals
P   <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]                                                            # The F0 individuals

onegenotypeF2 <- which(lapply(apply(genotypes[,F2], 1, table), length) < 2)                                               # Markers with only one genotype in the F2 cannot be used in QTL mapping
onegenotypeF1 <- which(lapply(apply(genotypes[,F1], 1, table), length) < 2)                                               # Markers with only one genotype in the F1
onegenotypeP  <- which(lapply(apply(genotypes[,P], 1, table), length) < 2)                                                # Markers with only one genotype in the P

genotypesF2   <- genotypes[-onegenotypeF2, F2]                                                                            # Only the unique markers in the F2
cat("Left with", nrow(genotypesF2), " F2 markers\n")                                                                      # == Left with 21260 markers    // Left with 11581 markers
genotypesF1   <- genotypes[-onegenotypeF1, F1]                                                                            # Only the unique markers in the F1
cat("Left with", nrow(genotypesF1), " F1 markers\n")                                                                      # == Left with 21260 markers    // Left with 11581 markers
genotypesP    <- genotypes[-onegenotypeP, P]                                                                              # Only the unique markers in the P
cat("Left with", nrow(genotypesP), " P markers\n")                                                                        # == Left with 21260 markers    // Left with 11581 markers

# Fixed effects
vater <- as.factor(phenotypes[F2, c("Vater")]) ; mutter <- as.factor(phenotypes[F2, c("Mutter")]) ; wsize <- as.numeric(phenotypes[F2, c("WG2")]) ; 
wlabel <- as.factor(phenotypes[F2, c("W.Label")]) ; season <- as.factor(phenotypes[F2, c("Season")]); fam <- as.factor(paste0(mutter, vater))

# Plot the different families
families <- lapply(unique(fam),function(x){ phenotypes[F2, "d63"][ which(fam == x) ] })
familiesG <- lapply(unique(fam),function(x){ genotypesF2["UNC5048297",][ which(fam == x) ] })
names(families) <- unique(fam)
names(familiesG) <- unique(fam)

F2toFam <- unlist(lapply(families, median))
sorder <- sort(F2toFam, index.return = TRUE)$ix
F2toFam <- F2toFam[ sorder ]

totalOffspring <- unlist(lapply(families[ names(F2toFam) ],length))

mothersF2 <- substr(names(F2toFam), 0, 8)                                 # Mothers of the F2
fathersF2 <- substr(names(F2toFam), 9, 16)                                # Dads of the F2

mGeno <- unlist(genotypesF1["UNC5048297", mothersF2])
fGeno <- unlist(genotypesF1["UNC5048297", fathersF2])
famGeno <- as.factor(paste0(mGeno, fGeno))

nA <- unlist(lapply(lapply(strsplit(as.character(famGeno),""), "==", "A"), sum))
nH <- unlist(lapply(lapply(strsplit(as.character(famGeno),""), "==", "H"), sum))
nB <- unlist(lapply(lapply(strsplit(as.character(famGeno),""), "==", "B"), sum))

plot(c(1, length(totalOffspring)),c(20, 65), t = 'n',xaxt='n', xlab="",ylab="Bodyweight (Grams)", main="Population structure", sub="Day 63")
boxplot(families[ names(F2toFam) ], col = rgb(1-(nB / 2), 1-(nA / 5), (nH / 2)), add=TRUE, las=2,cex.axis=0.2)
text(1:length(totalOffspring), rep(63,length(totalOffspring)), as.character(totalOffspring))
text(1:length(totalOffspring), 20, as.character(famGeno))

colN <- c("orange","gray","green")
names(colN) <- c("A","H","B")

for(x in 1:length(totalOffspring)){
  children <- rownames(phenotypes[F2,][fam==names(F2toFam)[x],])
  childG <- paste0(genotypesF2["UNC5048297", children],collapse="")
  gcolz <- colN[unlist(familiesG[ names(F2toFam) ][[x]])]
  points(rep(x, totalOffspring[x]), families[ names(F2toFam) ][[x]], col="black",pch=19, cex=1.1)
  points(rep(x, totalOffspring[x]), families[ names(F2toFam) ][[x]], col=gcolz,pch=18,cex=0.8)
  nA <- unlist(lapply(lapply(strsplit(as.character(childG),""), "==", "A"), sum))
  nH <- unlist(lapply(lapply(strsplit(as.character(childG),""), "==", "H"), sum))
  nB <- unlist(lapply(lapply(strsplit(as.character(childG),""), "==", "B"), sum))
  text(x-0.30, 61, paste0(nA),  cex=0.7, col="orange")
  text(x-0.15, 61, paste0("/"), cex=0.7)
  text(x+0.00, 61, paste0(nH),  cex=0.7)
  text(x+0.15, 61, paste0("/"), cex=0.7)
  text(x+0.30, 61, paste0(nB),  cex=0.7, col="green")
}




### Now do the adjustment of the phenotypes based on the BFMI locus in the children                                 
BFMIlocus <- as.numeric(as.factor(unlist(genotypesF2["UNC5048297", F2])))                                                   # Major effect locus
pheSubset <- phenotypes[F2, 5:36]                                                                                         # only the F2 individuals and the real phenotypes


# Fixed effect significance (shared factors)
significance <- t(apply(pheSubset, 2, function(x){ aa <- lm(as.numeric(x) ~ wlabel + season); return(anova(aa)[[5]]) }))[,-3]
colnames(significance) <- c("Litter number", "Season")
pheAdjusted <- apply(pheSubset, 2, function(x){ aa <- lm(as.numeric(x) ~ wlabel + season); 
                                         bb <- rep(NA, length(as.numeric(x))); 
                                         bb[as.numeric(names(aa$residuals))] <- (aa$residuals + aa$coefficients["(Intercept)"]); 
                                         return(bb) })

# Fixed effect significance (BFMI estimated factors)
sigBfmi <- t(apply(pheSubset[BFMIlocus == 1,], 2, function(x){ aa <- lm(as.numeric(x) ~ fam[BFMIlocus == 1] + wsize[BFMIlocus == 1]); return(anova(aa)[[5]]) }))[,-3]
colnames(sigBfmi) <- c("Family (BFMI)", "Litter size (BFMI)")
pheBfmi <- apply(pheSubset[BFMIlocus == 1,], 2, function(x){ aa <- lm(as.numeric(x) ~ fam[BFMIlocus == 1] + wsize[BFMIlocus == 1]);  
                                         bb <- rep(NA, length(as.numeric(x))); 
                                         bb[as.numeric(names(aa$residuals))] <- (aa$residuals + aa$coefficients["(Intercept)"]); 
                                         return(bb) })
rownames(pheBfmi) <- rownames(pheSubset[BFMIlocus == 1,])

# Fixed effect significance (B6+H estimated factors)
sigOthers <- t(apply(pheSubset[BFMIlocus != 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[BFMIlocus != 1] + wsize[BFMIlocus != 1]); return(anova(aa)[[5]]) }))[,-3]
colnames(sigOthers) <- c("Family (B6+H)", "Litter size (B6+H)")
pheOthers <- apply(pheSubset[BFMIlocus != 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[BFMIlocus != 1] + wsize[BFMIlocus != 1]); 
                                         bb <- rep(NA, length(as.numeric(x))); 
                                         bb[as.numeric(names(aa$residuals))] <- (aa$residuals + aa$coefficients["(Intercept)"]); 
                                         return(bb) })
rownames(pheOthers) <- rownames(pheSubset[BFMIlocus != 1,])

# Merge the different significances
fixedEffectSignificance <- cbind(significance, sigBfmi, sigOthers)

alldata           <- rbind(pheBfmi, pheOthers)              # Create the adjusted alldata
write.table(alldata, "Analysis/Phenotypes_newModelQTL.txt", sep = "\t",  quote=FALSE)
genotypesF2qtl    <- genotypesF2[,rownames(alldata)]        # Order has changed in the phenotypes, reorder the genotypes
write.table(genotypesF2qtl, "Analysis/Genotypes_newModelQTL.txt", sep = "\t",  quote=FALSE)

result <- matrix(NA, nrow(genotypesF2qtl), ncol(alldata))
for(x in 1:nrow(genotypesF2qtl)){
  for(y in 1:ncol(alldata)){
    tryCatch(res <- anova(lm(alldata[,y] ~ as.factor(as.character(genotypesF2qtl[x, ]))))[[5]], error = function(e){ res <<- rep(NA, 2) }) # Basic association test
    result[x, y] <- res[1]
  }
  cat("Done marker", x,"\n")
}
colnames(result) <- colnames(alldata)
rownames(result) <- rownames(genotypesF2qtl)
write.table(result, "Analysis/QTLresults_newModelQTL.txt", sep = "\t",  quote=FALSE)

for(x in 1:ncol(result)){
  plot(-log10(result[,x]), t = 'l', main=colnames(result)[x])
  scan()
}

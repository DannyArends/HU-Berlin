# QTL analysis of the QTL data from Mustapha
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

# TODO use forward selection

setwd("E:/Chicken/ClassicalPhenotypes/QTLanalysis_Mostafa_2010")

genotypes <- read.table("genotypes.txt", sep="\t", na.strings = ".", header=TRUE)                                 # Load the genotypes
bodycomposition <- read.table("bodycomposition.txt", sep="\t", na.strings = ".", header=TRUE)[, -1]               # Load the phenotypes

bodycomposition <- bodycomposition[]

BCinGT <- which(bodycomposition[,"ID"] %in% genotypes[,"ID"])                                                     # Which body composition traits have genotypes
bodycomposition <- bodycomposition[BCinGT, ]

GTinBC <- which(genotypes[,"ID"] %in% bodycomposition[,"ID"])                                                     # Which genotypes have body composition traits
genotypes <- genotypes[GTinBC, ]

GTorderToMatchBC <- match(as.character(bodycomposition[,"ID"]), as.character(genotypes[,"ID"]))                   # Reorder the genotypes such that they match the body composition ordering
genotypes <- genotypes[GTorderToMatchBC, ]

rownames(genotypes) <- genotypes[,"ID"]
rownames(bodycomposition) <- bodycomposition[,"ID"]

mapQTL <- function(bodycomposition, genotypes, phenotypename = "Body.weight.at.20.weeks..g"){                     # Function to do QTL mapping using 2 covariates
  markersForQTL <- colnames(genotypes)[7:ncol(genotypes)]                                                         # Markers that are used in QTL mapping
  LODscores <- NULL                                                                                               # Matrix for results
  for(markername in markersForQTL){
    m <- genotypes[,markername]
    m.dom <- as.numeric(m == 12)
    m.add <- (as.numeric(as.factor(m)) - 1)
    
    okM = which(!is.na(m.dom))
    bc <- bodycomposition[okM,]

    cof <- genotypes[okM,"UMA4.034"]
    
    m <- m[okM]
    m.dom <- m.dom[okM]
    m.add <- m.add[okM]
    
    mA <- lm(bc[,phenotypename] ~  bc[,"HatchDate"] + bc[,"Sire"] + bc[,"Cross"] + cof + as.factor(m))
    
    m1 <- lm(bc[,phenotypename] ~  bc[,"HatchDate"] + bc[,"Sire"] + bc[,"Cross"] + m.dom + m.add)
    m0 <- lm(bc[,phenotypename] ~  bc[,"HatchDate"] + bc[,"Sire"] + bc[,"Cross"])
    res  <- anova(m1)
    pval <- anova(m1, m0)
    ptest <- anova(mA, m0)
    LODscores <- rbind(LODscores, c(-log10(res[[5]]), -log10(pval[[6]][2]), -log10(ptest[[6]][2])))
  }
  rownames(LODscores) <- markersForQTL
  #colnames(LODscores) <- c("HatchDate", "Subfamily", "Cross", "M.Dom", "M.Add", "Residuals", "ModelDiff", "ModelDiff2")                                      # Add the row names to our result matrix
   return(LODscores)
}

phenotypesForQTL <- colnames(bodycomposition)[8:ncol(bodycomposition)]                                           # Phenotypes that are used in QTL mapping

LODresults <- vector("list", length(phenotypesForQTL))                                                           # Empty list of QTL matrices, one for each phenotype
names(LODresults) <- phenotypesForQTL
for(phenotype in phenotypesForQTL){
  LODresults[[phenotype]] <- mapQTL(bodycomposition, genotypes, phenotype)                                       # Map the QTL and store the resulting matrix in the list
}

plot(LODresults$"Body.weight.at.20.weeks..g"[,"ModelDiff"], t = 'l')
points(LODresults$"Body.weight.at.20.weeks..g"[,"M.Add"], t = 'l', col="green")
points(LODresults$"Body.weight.at.20.weeks..g"[,"M.Dom"], t = 'l', col='red')

for(x in 1:length(LODresults)){
  above4 <- which(LODresults[[x]][,"Marker"] > 4.65)                                                             # Analyse the resulting profiles, and look at which markers are above the threshold
  cat(names(LODresults)[x], names(above4),"\n")                                                                  # -log10(0.05/ (123 * 18)) = 4.65
}
plot(LODresults$Mesenteric.adipose.tissue..g[,"Marker"], t='l')                                                  # Plot one of our QTL profiles

## What do we learn: the first 6 traits show a similar QTL profile from:
## Small region: rs14488203 => 4 69583407 to rs14494169 => 4 78715886
## Extended region: MCW0276 => 4 59553432 to UMA4024 =>    4 84774762

topmarker <- "UMA4.034"


c1 <- which(bodycomposition[,"Cross"] == 1)   # NHI x WL77
c2 <- which(bodycomposition[,"Cross"] == 2)   # WL77 x NHI




ii <- which(!is.na(bc[,phenotypename]))
m <- genotypes[,topmarker]

pcor <- lm(bc[ii,phenotypename] ~  bc[ii,"HatchDate"] + bc[ii,"Sire"])$residuals

anova(lm(pcor ~ bc[ii,"Cross"] * as.factor(m[ii])))

anova(lm(bc[,phenotypename] ~  bc[,"HatchDate"] + bc[,"Sire"] + bc[,"Cross"] * as.factor(m)))



## Top marker "Body.weight.at.20.weeks..g" = rs14488203

phenotypename <- "Body.weight.at.20.weeks..g"
markername <- "rs14488203"

gt <- as.factor(as.character(genotypes[,markername]))
sire <- as.factor(as.character(bodycomposition[,"Sire"]))
hd <- as.factor(as.character(bodycomposition[,"HatchDate"]))


residual1 <- lm(bodycomposition[c1, "Body.weight.at.20.weeks..g"] ~ hd[c1] + sire[c1])$residuals
residual2 <- lm(bodycomposition[c2, "Body.weight.at.20.weeks..g"] ~ hd[c2] + sire[c2])$residuals


markersForQTL <- colnames(genotypes)[7:ncol(genotypes)]                                                         # Markers that are used in QTL mapping
LODscores <- NULL                                                                                               # Matrix for results
for(markername in markersForQTL){
  res <- anova(lm(residual ~  bodycomposition[as.numeric(names(residual)),"Cross"] * genotypes[as.numeric(names(residual)),markername]))
  LODscores <- rbind(LODscores, -log10(res[[5]]))
}
rownames(LODscores) <- markersForQTL


markername <- "UMA4.034"


ma <- anova(lm(bodycomposition[,phenotypename] ~  hd + sire + gt))
ma[[2]] / sum(ma[[2]])

AIC(lm(bodycomposition[,phenotypename] ~ gt))               # Best model: 3782.041
AIC(lm(bodycomposition[,phenotypename] ~ hd + gt))          # better model: 3754.228
AIC(lm(bodycomposition[,phenotypename] ~ hd + sire + gt))   # better model: 3736.432

residual <- rep(NA, length(bodycomposition[,phenotypename]))
names(residual) <- 1:length(residual)
ressss <- lm(bodycomposition[,phenotypename] ~ hd + sire)$residuals
residual[names(ressss)] <- ressss
plot(residual ~ gt)



m <- genotypes[,markername]
m.dom <- as.numeric(m == 12)
m.add <- (as.numeric(as.factor(m)) - 1)

anova(lm(bodycomposition[,phenotypename] ~  bodycomposition[,"HatchDate"] + bodycomposition[,"Sire"] + m))
anova(lm(bodycomposition[,phenotypename] ~  bodycomposition[,"HatchDate"] + bodycomposition[,"Sire"] + m.add))


okM = which(!is.na(m.dom))
bodycomposition <- bodycomposition[okM,]

m.dom <- m.dom[okM]
m.add <- m.add[okM]

model1 <- lm(bodycomposition[,phenotypename] ~  bodycomposition[,"HatchDate"] + bodycomposition[,"Sire"] + m.dom + m.add)
model0 <- lm(bodycomposition[,phenotypename] ~  bodycomposition[,"HatchDate"] + bodycomposition[,"Sire"])
anova(model1,model0)
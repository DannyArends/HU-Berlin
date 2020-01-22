#
# Analyze the butana population for Somatic Cell Score
#

getSeason <- function(mmonths) {
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 8] <- "DRY"
  ret[mmonths >= 9 & mmonths <= 12 | mmonths == 1 | mmonths == 2] <- "WET"
  return(ret)
}

library(lme4)

setwd("D:/Edrive/Cow/Salma Butana")
butanadata <- read.csv("phenotypes_2020.txt", sep = "\t", header=TRUE, check.names=FALSE)
butana <- butanadata[which(butanadata[, "Breed"] == "Butana"),]

# Get the phenotype we are analyzing
phenotype <- log2(as.numeric(as.character(butana[, "SCC"])) / 100000) + 3

# Get the covariates we want to investigate in model building
animal <- butana[, "ID"]
lactation <- as.factor(butana[, "Number of calvings"] - 1)
year <- unlist(lapply(strsplit(as.character(butana[,"Sampledate"]), "/"), "[", 3))
birthyear <- unlist(lapply(strsplit(as.character(butana[,"Birthdate"]), "/"), "[", 3))
month <- as.numeric(unlist(lapply(strsplit(as.character(butana[,"Sampledate"]), "/"), "[", 1)))
season <- getSeason(month)
daysinmilk <- butana[, "Lactation length"]
firstCalf <- butana[, "Age at first calving"]

# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(phenotype), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Year = as.factor(year),
                    Birthyear = as.factor(birthyear),
                    Season = as.factor(season),
                    Firstcalf = as.factor(firstCalf), 
                    DIM = as.numeric(daysinmilk))

# Remove rows with missing data
hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }

# Which covariates influence the phenotype in Butana ?
null.model <- lmer(Y ~ (1|Animal), data = data.frame(mdata), REML = FALSE)
model1 <- lmer(Y ~ Lact + (1|Animal), data = data.frame(mdata), REML = FALSE)
anova(model1, null.model) # N.S. p = 0.2931

model2 <- lmer(Y ~ Season + (1|Animal), data = data.frame(mdata), REML = FALSE)
anova(model2, null.model) # N.S. p = 0.506

model3 <- lmer(Y ~ Year + (1|Animal), data = data.frame(mdata), REML = FALSE)
anova(model3, null.model) # S p = 0.02678

model4 <- lmer(Y ~ Birthyear + (1|Animal), data = data.frame(mdata), REML = FALSE)
anova(model4, null.model) # S p = 0.02842

model5 <- lmer(Y ~ Firstcalf + (1|Animal), data = data.frame(mdata), REML = FALSE)
anova(model5, null.model) # S p = 2.505e-05

model6 <- lmer(Y ~ DIM + (1|Animal), data = data.frame(mdata), REML = FALSE)
anova(model6, null.model) # N.S. p = 0.9088

# Get the marker data as genotypes and perform some QC on it
genotypes <- butana[,29:ncol(butana)]
for (x in 1:ncol(genotypes)) {
  tbl <- table(genotypes[,x])
  for (gt in names(tbl)) {
    nAnimals <- length(unique(butana[which(genotypes[,x] == gt), "ID"]))
    if (nAnimals < 10) {
      genotypes[which(genotypes[,x] == gt), x] <- NA
    }
  }
}
seggregates <- which(unlist(lapply(apply(genotypes,2,table), length)) > 1)
genotypes <- genotypes[, seggregates]

pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata <- cbind(Y = phenotype, 
                 Marker = genotypes[,x],
                 Birthyear = as.factor(birthyear),
                 Year = as.factor(year), 
                 Firstcalf = as.factor(firstCalf),                  
                 Animal = animal)

  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- lmer(Y ~ Firstcalf + Year + Birthyear + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Firstcalf + Year + Birthyear + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
}
names(pvals) <- colnames(genotypes)
pvals

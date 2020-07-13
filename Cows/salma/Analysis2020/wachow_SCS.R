#
# Analyze the holstein population for Wachow for Somatic Cell Score
#

getSeason <- function(mmonths) {
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5] <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8] <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11] <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2] <- "Winter"
  return(ret)
}

library(lme4)

setwd("D:/Edrive/Cow/Salma Wachow")
wachowdata <- read.csv("phenotypes2020.wachow_v2.txt", sep = "\t", header=TRUE, check.names=FALSE)

# Get the phenotype we are analyzing
phenotype <- log2(as.numeric(as.character(wachowdata[, "SCC"])) / 100) + 3
phenotype[which(!is.finite(phenotype))] <- 0

# Get the covariates we want to investigate in model building
animal <- wachowdata[, "LOM"]
lactation <- as.factor(wachowdata[, "Lactation"])
year <- unlist(lapply(strsplit(as.character(wachowdata[,"Sampledate"]), "-"), "[", 1))
birthyear <- unlist(lapply(strsplit(as.character(wachowdata[,"Birthdate"]), "-"), "[", 1))
month <- as.numeric(unlist(lapply(strsplit(as.character(wachowdata[,"Sampledate"]), "-"), "[", 2)))
season <- getSeason(month)
firstCalf <- wachowdata[, "FirstCalfIndays"]

# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(phenotype), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Year = as.factor(year),
                    Birthyear = as.factor(birthyear),
                    Season = as.factor(season),
                    Firstcalf = as.factor(firstCalf))

# Remove rows with missing data
hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }

# Which covariates influence the phenotype in Butana ?
null.model <- lmer(Y ~ (1|Animal), data = mdata, REML = FALSE)
model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata, REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Season + (1|Animal), data = mdata, REML = FALSE)
anova(model2, null.model) # 4.425e-07 ***

model3 <- lmer(Y ~ Year + (1|Animal), data = mdata, REML = FALSE)
anova(model3, null.model) # < 2.2e-16 ***

model4 <- lmer(Y ~ Birthyear + (1|Animal), data = mdata, REML = FALSE)
anova(model4, null.model) # 3.12e-05 ***

model5 <- lmer(Y ~ Firstcalf + (1|Animal), data = mdata, REML = FALSE)
anova(model5, null.model) # 0.3098

genotypes <- wachowdata[,10:ncol(wachowdata)]
for (x in 1:ncol(genotypes)) {
  tbl <- table(genotypes[,x])
  for (gt in names(tbl)) {
    nAnimals <- length(unique(wachowdata[which(genotypes[,x] == gt), "ID"]))
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
                 Lact = as.factor(lactation),
                 Season = as.factor(season),                 
                 Animal = animal)

  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- lmer(Y ~ Year + Lact + Birthyear + Season + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Year + Lact + Birthyear + Season + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
}
names(pvals) <- colnames(genotypes)
pvals

#R.HAS.M4 = suggestive = 0.068160937
#R.HAS.M10 = significant = 0.002713075
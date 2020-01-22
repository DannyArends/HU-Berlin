#
# Analyze the cross breed population for Somatic Cell Score
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
crossbreed <- butanadata[which(butanadata[, "Breed"] != "Butana"),]

crossbreed[which(crossbreed[, "SCC"] == 0), "SCC"] <- NA

# Get the phenotype we are analyzing
phenotype <- log2(as.numeric(as.character(crossbreed[, "SCC"])) / 100000) + 3

# Get the covariates we want to investigate in model building
animal <- crossbreed[, "ID"]
lactation <- as.factor(crossbreed[, "Number of calvings"] - 1)
farm <- crossbreed[, "Farm"]
year <- unlist(lapply(strsplit(as.character(crossbreed[,"Sampledate"]), "/"), "[", 3))
birthyear <- unlist(lapply(strsplit(as.character(crossbreed[,"Birthdate"]), "/"), "[", 3))
month <- as.numeric(unlist(lapply(strsplit(as.character(crossbreed[,"Sampledate"]), "/"), "[", 1)))
season <- getSeason(month)
daysinmilk <- crossbreed[, "Lactation length"]
firstCalf <- crossbreed[, "Age at first calving"]

# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(phenotype), 
                    Lact = as.factor(lactation), 
                    Animal = as.factor(animal), 
                    Farm = as.factor(farm), 
                    Year = as.factor(year), 
                    Birthyear = as.factor(birthyear), 
                    Season = as.factor(season), 
                    Firstcalf = as.factor(firstCalf), 
                    DIM = as.numeric(daysinmilk))
                    
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
if(length(hasMissing) > 0){ mdata <- mdata[-hasMissing,] }

### Model different covariates
null.model <- lmer(Y ~ (1|Animal), data = mdata, REML = FALSE)
m1 <- lmer(Y ~ Lact + (1|Animal), data = mdata, REML = FALSE)
anova(m1, null.model) # N.S. p = 0.4919

m2 <- lmer(Y ~ Farm + (1|Animal), data = mdata, REML = FALSE)
anova(m2, null.model) # Sug. p = 0.06439 (include anyway)

m3 <- lmer(Y ~ Year + (1|Animal), data = mdata, REML = FALSE)
anova(m3, null.model) # N.S. p = 0.107

m4 <- lmer(Y ~ Birthyear + (1|Animal), data = mdata, REML = FALSE)
anova(m4, null.model) # N.S. p = 0.2008

m5 <- lmer(Y ~ Season + (1|Animal), data = mdata, REML = FALSE)
anova(m5, null.model) # N.S. p = 0.7446

m6 <- lmer(Y ~ Firstcalf + (1|Animal), data = mdata, REML = FALSE)
anova(m6, null.model) # S. p = 0.02675

m7 <- lmer(Y ~ DIM + (1|Animal), data = mdata, REML = FALSE)
anova(m7, null.model) # N.S. p = 0.1884

# Get the marker data as genotypes and perform some QC on it
genotypes <- crossbreed[,29:ncol(crossbreed)]
for (x in 1:ncol(genotypes)) {
  tbl <- table(genotypes[,x])
  for(gt in names(tbl)){
    nAnimals <- length(unique(crossbreed[which(genotypes[,x] == gt), "ID"]))
    if(nAnimals < 10){
      genotypes[which(genotypes[,x] == gt), x] <- NA
    }
  }
}
seggregates <- which(unlist(lapply(apply(genotypes,2,table), length)) > 1)
genotypes <- genotypes[, seggregates]

pvals <- NULL
for(x in 1:ncol(genotypes)) {
  # Recombine all data into a data.frame add the genotypes
  mdata <- data.frame(Y = as.numeric(phenotype), 
                      Marker = as.factor(genotypes[,x]),
                      Animal = animal,
                      Farm = as.factor(farm),
                      Firstcalf = as.factor(firstCalf))

  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if(length(hasMissing) > 0){ mdata <- mdata[-hasMissing,] }
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- lmer(Y ~ Farm + Firstcalf + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Farm + Firstcalf + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
}
names(pvals) <- colnames(genotypes)
pvals

# Analysis of Sudanese Goat Breeds
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014

setwd("D:/Siham/")
library(lme4)
goatdata <- read.table("Sudan_Goat_Breeds.txt",sep="\t", header=TRUE,na.strings="")

genotypes <- c("DGAT1_K232A", "DGAT1_VNTR", "Leptin", "K_casein")
traits <- c("MY", "Fat", "Protein", "Fat_kg", "Protein_Kg")

goatdata <- goatdata[goatdata[,"Breed"]==1,]

goatdata <- goatdata[-which(is.na(goatdata[,"Year_of_calving"])),]

resulttable <- NULL
for(trait in traits){
  for(genotype in genotypes){
    idxNA <- which(is.na(goatdata[,genotype]))
    if(length(idxNA) > 1){
      gdata <- goatdata[-idxNA,]
    }else{
      gdata <- goatdata
    }
    traitvalues <- gdata[,trait]
    genotypevalues <- as.numeric(gdata[,genotype])
    X1 <- (gdata[,"DIM"]/305)
    X2 <- X1*X1
    X3 <- log(305 / gdata[,"DIM"])
    X4 <- X3*X3
    
    model <- lmer(traitvalues ~ genotypevalues + (1 | gdata[,"Proben_No"]) )
    coefs <- data.frame(coef(summary(model)))
    pval <- 2 * (1 - pnorm(abs(coefs$t.value)))
    resulttable <- rbind(resulttable, c(trait, genotype, pval))
  }
}
#colnames(resulttable) <- c("Trait", "Marker", "P-value (raw)")
resulttable <- as.data.frame(resulttable)
resulttable

#write.table(resulttable,file="AnalysisResultsB1.txt",sep="\t",row.names=FALSE)
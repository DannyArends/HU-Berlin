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

#goatdata <- goatdata[goatdata[,"Breed"]==3,]

resulttable <- NULL
for(trait in traits){
  for(genotype in genotypes){
    avgtraitvalue <- NULL
    avgDIM <- NULL
    for(x in unique(goatdata[,"Proben_No"])){
      idxs <- which(goatdata[,"Proben_No"] == x)
      avgtraitvalue <- c(avgtraitvalue, mean(goatdata[idxs, trait], na.rm=TRUE))
      avgDIM <- c(avgDIM, mean(goatdata[idxs, "DIM"], na.rm=TRUE))
    }
  
    gdata <- goatdata[!duplicated(goatdata[,"Proben_No"]),]

    genotypevalues <- as.numeric(gdata[,genotype])
    
    model1 <- lm(avgtraitvalue ~ gdata[,"Breed"] + gdata[,"Year_of_calving"] + gdata[,"Season_of_calving"] + avgDIM + gdata[,"Lact"] + genotypevalues )
    model2 <- lm(avgtraitvalue ~ gdata[,"Breed"] + genotypevalues )
    
    resulttable <- rbind(resulttable, c(trait, genotype, anova(model2)[[5]][2], anova(model1)[[5]][6], as.numeric(anova(model1)[[5]][6] < 0.0025)))
  }
}
colnames(resulttable) <- c("Trait", "Marker", "P-value (raw)", "P-value (cor)", "isSignificant")
resulttable <- as.data.frame(resulttable)
resulttable

#write.table(resulttable,file="AnalysisResults_ALL.txt",sep="\t",row.names=FALSE)




# RNA Seq - Classical phenotypes analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2014
# first written Nov, 2014

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")
F1phe  <- read.csv("20140718_F1V.txt", sep="\t", header=TRUE, check.names=FALSE)

F1pheFF <- F1phe[which(F1phe[,"futter"] == "FF"),]

F1pheFF_M <- F1pheFF[which(F1pheFF[,"sex"] == "m"),]
F1pheFF_F <- F1pheFF[which(F1pheFF[,"sex"] == "f"),]

# d21   d28   d35   d42   d49   d56   d63   d70
cat("day\tbw\tfat\tlean\tfat/lean\n", file="Analysis/statistics_males.txt")
for(x in c("d28","d35","d42","d49","d56","d63","d70")){
  bwa <- anova(lm(F1pheFF_M[,x]                                                   ~ F1pheFF_M[,"WG"] + F1pheFF_M[,"dirocross"]))
  fat <- anova(lm(F1pheFF_M[,paste0(x,"mrifat")]                                  ~ F1pheFF_M[,"WG"] + F1pheFF_M[,"dirocross"]))
  lean <- anova(lm(F1pheFF_M[,paste0(x,"mrilean")]                                ~ F1pheFF_M[,"WG"] + F1pheFF_M[,"dirocross"]))
  fl <- anova(lm(F1pheFF_M[,paste0(x,"mrifat")] / F1pheFF_M[,paste0(x,"mrilean")] ~ F1pheFF_M[,"WG"] + F1pheFF_M[,"dirocross"]))
  cat(x, bwa[[5]][2], fat[[5]][2], lean[[5]][2], fl[[5]][2], "\n", sep="\t", file="Analysis/statistics_males.txt", append=TRUE)
}

# d21   d28   d35   d42   d49   d56   d63   d70
cat("day\tbw\tfat\tlean\tfat/lean\n", file="Analysis/statistics_females.txt")
for(x in c("d28","d35","d42","d49","d56","d63","d70")){
  bwa <- anova(lm(F1pheFF_F[,x]                                                   ~ F1pheFF_F[,"WG"] + F1pheFF_F[,"dirocross"]))
  fat <- anova(lm(F1pheFF_F[,paste0(x,"mrifat")]                                  ~ F1pheFF_F[,"WG"] + F1pheFF_F[,"dirocross"]))
  lean <- anova(lm(F1pheFF_F[,paste0(x,"mrilean")]                                ~ F1pheFF_F[,"WG"] + F1pheFF_F[,"dirocross"]))
  fl <- anova(lm(F1pheFF_F[,paste0(x,"mrifat")] / F1pheFF_F[,paste0(x,"mrilean")] ~ F1pheFF_F[,"WG"] + F1pheFF_F[,"dirocross"]))
  cat(x, bwa[[5]][2], fat[[5]][2], lean[[5]][2], fl[[5]][2], "\n", sep="\t", file="Analysis/statistics_females.txt", append=TRUE)
}

cat("phenotype\tmale\tfemale\n", file="Analysis/statistics_invasive.txt")
for(ph in c("GF1", "GF2", "totalGF", "RF1", "RF2", "totalRF", "IF", "muskel", "leber", "BAT", "LD")){
  male   <- anova(lm(F1pheFF_M[, ph] / F1pheFF_M[, "d70"]  ~ F1pheFF_M[,"WG"] + F1pheFF_M[,"dirocross"]))
  female <- anova(lm(F1pheFF_F[, ph] / F1pheFF_F[, "d70"]  ~ F1pheFF_F[,"WG"] + F1pheFF_F[,"dirocross"]))
  cat(paste0(ph,"/BW(70)"), male[[5]][2], female[[5]][2],"\n",sep="\t", file="Analysis/statistics_invasive.txt", append=TRUE)
}


boxplot(F1pheFF_M[F1pheFF_M[,"dirocross"] == "matBFMI",c("d28","d35","d42","d49","d56","d63","d70")], col="orange")
boxplot(F1pheFF_M[F1pheFF_M[,"dirocross"] == "matB6",c("d28","d35","d42","d49","d56","d63","d70")], add=TRUE, col="gray")

boxplot(F1pheFF_M[F1pheFF_M[,"dirocross"] == "matBFMI",paste0(c("d28","d35","d42","d49","d56","d63","d70"),"mrifat")], col="orange")
boxplot(F1pheFF_M[F1pheFF_M[,"dirocross"] == "matB6",paste0(c("d28","d35","d42","d49","d56","d63","d70"),"mrifat")], add=TRUE, col="gray")

 

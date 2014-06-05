# analysisOfClassicalPhenotypesInHorses.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Analysis to create some plots for Saria's CTC 2014 poster

setwd("D:/Saria_Horses")
rawdata <- read.table("A.Pferde Projekt.Saria.csv", sep=';', header=TRUE)
interesting <- c("Saklawi","Kahlawi","Hamdani")

#Subset, take the horses were intrested in
rows <- which(rawdata[,"Rassan"] %in% interesting)
rawdata <- rawdata[rows,]

#Ages 2 groups
ages <- rawdata[,"D.Birth"]
agegroup <- as.numeric(ages < 2007)

#Race, sex and phenotype names
rassen <- as.character(rawdata[,"Rassan"])
sex <- rawdata[,"Sex"]
phenotypes <- c("WH","CW","CH","NG","TG","ChG","ChD","ChW","BLL","BL","FCL","HCL")

#Subset, take the phenotypes were intrested in
phenotypedata <- rawdata[,phenotypes]

#Some basic statistics
mm <- NULL
dd <- NULL
si <- NULL
for(phe in phenotypes){
  result <- lm(phenotypedata[,phe] ~ as.factor(rassen) + agegroup + sex)
  si  <- rbind(si, anova(result)[[5]][2:3])
  plot(result$residuals ~ as.factor(rassen), main=phe)
  column <- NULL
  columndd <- NULL
  # Histogram of all phenotype data
  jpeg(paste0("Histograms/hist_", phe, ".jpg"))
    hist(phenotypedata[,phe], main = phe)
  dev.off()
  for(ras in interesting){
    ttest <- t.test(result$residuals[rassen==ras], result$residuals)
    column <- c(column, -log10(ttest$p.value))
    columndd <- c(columndd, as.numeric(ttest$estimate[1] > ttest$estimate[2]))
    
    #Histogram phenotype data per race
    jpeg(paste0("Histograms/hist_",ras,"_", phe, ".jpg"))
      hist(phenotypedata[which(rassen==ras), phe], main = paste0(phe," ",ras))
    dev.off()
  }
  mm <- rbind(mm, column)
  dd <- rbind(dd, columndd)
}

#Name the columns and rows
colnames(mm) <- interesting
rownames(mm) <- phenotypes

dd[dd==0] <- -1

colnames(dd) <- interesting
rownames(dd) <- phenotypes

colnames(si) <- c("agegroup","sex")
rownames(si) <- phenotypes
si <- -log10(si)

library(reshape2)
library(ggplot2)
library(plyr)
library(scales)

# Plot the data per race
mm.m <- melt(mm)
colnames(mm.m) <- c("Phenotype", "Rassan", "LOD")
mm.m <- ddply(mm.m, .(Phenotype), transform) #, rescale = rescale(value))
jpeg(file="heatmap.jpg", width = 1280, height = 786)
(p <- ggplot(mm.m, aes(Phenotype, Rassan)) + geom_tile(aes(fill = LOD), colour = "white") + scale_fill_gradient(low = "white", high = "red"))
p + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))+ theme(axis.ticks = element_blank()) #legend.position = "none"
dev.off()

# Plot the data per races, scaled for the direction
mm.m <- melt(mm * dd)
colnames(mm.m) <- c("Phenotype", "Rassan", "LOD")
mm.m <- ddply(mm.m, .(Phenotype), transform) #, rescale = rescale(value))
jpeg(file="heatmap_direction.jpg", width = 1280, height = 786)
(p <- ggplot(mm.m, aes(Phenotype, Rassan)) + geom_tile(aes(fill = LOD), colour = "white") + scale_fill_gradient(low = "yellow", high = "violet"))
p + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))+ theme(axis.ticks = element_blank()) #legend.position = "none"
dev.off()

#Plot the class data per phenotype
si.m <- melt(si)
colnames(si.m) <- c("Phenotype", "Class", "LOD")
si.m <- ddply(si.m, .(Phenotype), transform) #, rescale = rescale(value))
jpeg(file="heatmap_classes.jpg", width = 1280, height = 786)
(p <- ggplot(si.m, aes(Phenotype, Class)) + geom_tile(aes(fill = LOD), colour = "white") + scale_fill_gradient(low = "white", high = "red"))
p + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))+ theme(axis.ticks = element_blank()) #legend.position = "none"
dev.off()
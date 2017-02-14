# analysis.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2014
# first written Nov, 2014
#

install.packages("extrafont")
library("extrafont")
font_import()
fonts()

library("extrafont")
setwd("D:/Ddrive/Collegues/Ramona/")

rawdata <- read.table("Period_FA_EU.txt", sep="\t", header=TRUE)

anova(lm(rawdata[,"C12"] ~ as.factor(rawdata[,"Period"]) + rawdata[,"Group_EU"]))

anova(lm(rawdata[which(rawdata[,"Period"] == 1),"C12"] ~  rawdata[which(rawdata[,"Period"] == 1),"Group_EU"]))
anova(lm(rawdata[which(rawdata[,"Period"] == 2),"C12"] ~  rawdata[which(rawdata[,"Period"] == 2),"Group_EU"]))
anova(lm(rawdata[which(rawdata[,"Period"] == 3),"C12"] ~  rawdata[which(rawdata[,"Period"] == 3),"Group_EU"]))

box <- read.table("data_test.txt", sep = "\t", header=TRUE)
pval <- read.table("dana.pvalues.txt", sep = "\t", header=TRUE,row.names=1)


names(box)

plotMe <- function(box, pval, what="C12", main="MAIN"){
  # Create groups
  lowCalf6 <- box[box[,"G_comp"] == 1 & box[,"period"] == 1,]
  lowCalf12 <- box[box[,"G_comp"] == 1 & box[,"period"] == 2,]
  highCalf6 <- box[box[,"G_comp"] == 2 & box[,"period"] == 1,]
  highCalf12 <- box[box[,"G_comp"] == 2 & box[,"period"] == 2,]
  cow <- box[box[,"G_comp"] == 3 & box[,"period"] == 4,]

  # Create the plot
  plot(x=c(0.7,5.5),y = c(0,1.4 * max(box[,what])),ylab=main,t ='n',xaxt='n',xlab="Group / Feed efficiency", main=main)
  rect(-0.5,-0.5, 2.55, 1.4 * max(box[,what]) + 2.5, col=rgb(0.9,0.9,0.9,0.5), border = NA)
  rect(2.55,-0.5, 4.45, 1.4 * max(box[,what]) + 2.5, col=rgb(0.5,0.5,0.5,0.5), border = NA)
  boxplot(lowCalf6[,what], at=1,add=TRUE, col=rgb(0.8,0.8,0.8,1))
  boxplot(highCalf6[,what], at=2,add=TRUE, col=rgb(0.2,0.2,0.2,0.5))
  boxplot(lowCalf12[,what], at=3,add=TRUE, col=rgb(0.8,0.8,0.8,1))
  boxplot(highCalf12[,what], at=4,add=TRUE, col=rgb(0.2,0.2,0.2,0.5))
  boxplot(cow[,what],at=5,add=TRUE)
  axis(1,at=1:5, c("Low FE", "High FE", "Low FE", "High FE", "Cow"))
  text(1.5,1.4 * max(box[,what]),"Calf\nWeek 6")
  text(3.5,1.4 * max(box[,what]),"Calf\nWeek 12")
  text(5,1.4 * max(box[,what]),"Cow\nWeek 6")
  
  threshold <- c((0.05/31), (0.01/31), (0.0001/31))
  height <- 1.1 * max(box[,what])
  if(pval[what, 1] < threshold){  lines(x=c(1,2), y = c(height, height), lty=1); text(1.5,height + 0.1,c("*","**","***")[sum(pval[what,1] < threshold)],cex=1.1); height <- height + 0.25 }
  if(pval[what, 2] < threshold){  lines(x=c(1,3), y = c(height, height), lty=1); text(2,height + 0.1,c("*","**","***")[sum(pval[what,2] < threshold)],cex=1.1); height <- height + 0.25 }
  if(pval[what, 3] < threshold){  lines(x=c(1,4), y = c(height, height), lty=1); text(2.5,height + 0.1,c("*","**","***")[sum(pval[what,3] < threshold)],cex=1.1); height <- height + 0.25 }
  if(pval[what, 4] < threshold){  lines(x=c(1,5), y = c(height, height), lty=1); text(3,height + 0.1,c("*","**","***")[sum(pval[what,4] < threshold)],cex=1.1); height <- height + 0.25 }
  
  if(pval[what, 5] < threshold){  lines(x=c(2,3), y = c(height, height), lty=1); text(2.5,height + 0.1,c("*","**","***")[sum(pval[what,5] < threshold)],cex=1.1); height <- height + 0.25 }
  if(pval[what, 6] < threshold){  lines(x=c(2,4), y = c(height, height), lty=1); text(3,height + 0.1,c("*","**","***")[sum(pval[what,6] < threshold)],cex=1.1); height <- height + 0.25 }
  if(pval[what, 7] < threshold){  lines(x=c(2,5), y = c(height, height), lty=1); text(3.5,height + 0.1,c("*","**","***")[sum(pval[what,7] < threshold)],cex=1.1); height <- height + 0.25 }

  if(pval[what, 8] < threshold){  lines(x=c(3,4), y = c(height, height), lty=1); text(3.5,height + 0.1,c("*","**","***")[sum(pval[what,8] < threshold)],cex=1.1); height <- height + 0.25 }
  if(pval[what, 9] < threshold){  lines(x=c(3,5), y = c(height, height), lty=1); text(4,height + 0.1,c("*","**","***")[sum(pval[what,9] < threshold)],cex=1.1); height <- height + 0.25 }
  
  if(pval[what, 10] < threshold){ lines(x=c(4,5), y = c(height, height), lty=1); text(4.5,height + 0.1,c("*","**","***")[sum(pval[what,10] < threshold)],cex=1.1); height <- height + 0.25 }
}

plotMe(box, pval,"C12", "C12:0 content (%) in hair")

plotMe(box, pval,"C18D1t11", "C18:1trans-11 content (%) in hair")

plotMe(box, pval,"C18D2cis", "C18:2n-6 content (%) in hair")


# analysis.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2014
# first written Nov, 2014
#

setwd("d:/Ramona")

rawdata <- read.table("Period_FA_EU.txt", sep="\t", header=TRUE)

anova(lm(rawdata[,"C12"] ~ as.factor(rawdata[,"Period"]) + rawdata[,"Group_EU"]))

anova(lm(rawdata[which(rawdata[,"Period"] == 1),"C12"] ~  rawdata[which(rawdata[,"Period"] == 1),"Group_EU"]))
anova(lm(rawdata[which(rawdata[,"Period"] == 2),"C12"] ~  rawdata[which(rawdata[,"Period"] == 2),"Group_EU"]))
anova(lm(rawdata[which(rawdata[,"Period"] == 3),"C12"] ~  rawdata[which(rawdata[,"Period"] == 3),"Group_EU"]))
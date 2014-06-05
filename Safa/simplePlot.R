# simplePlot.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Basic line plot, for cattle data

setwd("D:/safa")
cattlenumbers <- read.csv("cattlenumber.csv")

plot(cattlenumbers, type = 'l', col = "red", lwd = 3)


KandBcattlenumber <- read.table("KandBcattlenumber.txt", row.names=1)
plot(c(1,2),c(0,1100000), type="n", col="blue", lwd=2,xaxt='n',xlab="Year", ylab="Cattle Number", yaxt='n')
points(KandBcattlenumber[,1], t='l', lwd=3)
points(KandBcattlenumber[,2], t='l', lwd=3, col="red")
axis(1, at=c(1,2), c(1999,2010))
axis(2, at=c(0,250000,500000,750000,1000000), c("0","250.000","500.000","750.000","1000.000"))
legend("bottomleft", colnames(KandBcattlenumber), col=c("black","red"), lwd=2)

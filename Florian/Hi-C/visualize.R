setwd("C:/Users/Arends/Desktop")

mBFMI <- read.table("subset1.binned.bfmi.alignments.txt")
mB6N <- read.table("subset1.binned.b6n.alignments.txt")

#Creating new chromosome 1 at bin: 0
#Creating new chromosome 10 at bin: 480062
#Creating new chromosome 11 at bin: 794518
#Creating new chromosome 12 at bin: 1089974
#Creating new chromosome 13 at bin: 1379283
#Creating new chromosome 14 at bin: 1671486
#Creating new chromosome 15 at bin: 1974394
#Creating new chromosome 16 at bin: 2223894
#Creating new chromosome 17 at bin: 2456144
#Creating new chromosome 18 at bin: 2685237
#Creating new chromosome 19 at bin: 2900852
#Creating new chromosome 2 at bin: 3044277
#Creating new chromosome 3 at bin: 3485695
#Creating new chromosome 4 at bin: 3877072
#Creating new chromosome 5 at bin: 4261738
#Creating new chromosome 6 at bin: 4630138
#Creating new chromosome 7 at bin: 4994893
#Creating new chromosome 8 at bin: 5365530
#Creating new chromosome 9 at bin: 5676520
#Creating new chromosome MT at bin: 5978627
#Creating new chromosome X at bin: 5978663
#Creating new chromosome Y at bin: 6417860 (6669677)


sChr3 <- 3485695
eChr3 <- 3877072

iix1 <- c(which(mBFMI[,2] > sChr3 & mBFMI[,2] < eChr3), which(mBFMI[,3] > sChr3 & mBFMI[,3] < eChr3))
iix2 <- c(which(mB6N[,2] > sChr3 & mB6N[,2] < eChr3), which(mB6N[,3] > sChr3 & mB6N[,3] < eChr3))

mBFMI <- mBFMI[iix1, ]
mB6N <- mB6N[iix2, ]

#plot(c(sChr3,eChr3), c(0, 6669677), t = 'n')
#op <- par(mfrow = c(2,1))
plot( c(0, 6669677), c(sChr3,eChr3), t = 'n', main = "BFMI Chr3")
points(mBFMI[,2], mBFMI[,3], pch=19, cex=0.1, col = rgb(0,0,0,0.01))

plot(c(sChr3,eChr3), c(sChr3,eChr3), t = 'n', main = "B6N Chr3")
points(mB6N[,2], mB6N[,3], pch=19, cex=0.1, col = rgb(0,0,0,0.01))


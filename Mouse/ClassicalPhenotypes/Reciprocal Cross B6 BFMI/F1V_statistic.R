# F1V_statistic.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Stefan Kärst, Danny Arends
# last modified Aug, 2014
# first written Jul, 2014

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")
f <- read.table(file="20140718_F1V.txt", sep="\t", na.strings = "NA", stringsAsFactors=F, header=TRUE )

table(f$futter == "FF" & f$sex =="f" & f$dirocross=="matB6")
table(f$futter == "ZF" & f$sex =="f")

summ <- summary(f[,1:78])

table(unique(mb6$mutter))
table(unique(mb6$vater))

head(f)

mb6  <- subset(f, f$dirocross == "matB6"   & f$futter == "FF");                 cat("Number of animals (m+f, maternalB6, Fat Food):", nrow(mb6), "\n")
mbf  <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF");                 cat("Number of animals (m+f, maternalBFMI, Fat Food):", nrow(mbf), "\n")
mb6m <- subset(f, f$dirocross == "matB6"   & f$futter == "FF" & f$sex== "m");   cat("Number of animals (male, maternalB6, Fat Food):", nrow(mb6m), "\n")
mbfm <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF" & f$sex== "m");   cat("Number of animals (male, maternalBFMI, Fat Food):", nrow(mbfm), "\n")
mb6f <- subset(f, f$dirocross == "matB6"   & f$futter == "FF" & f$sex== "f");   cat("Number of animals (female, maternalB6, Fat Food):", nrow(mb6f), "\n")
mbff <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF" & f$sex== "f");   cat("Number of animals (female, maternalBFMI, Fat Food):", nrow(mbff), "\n")

summb6 <- summary(mb6[,1:78])
summbf <- summary(mbf[,1:78])

## test auf normalverteilung
shapiro.test(mb6$d21)       # n.nd
shapiro.test(mb6m$d21)      # nd
shapiro.test(mb6f$d21)      # n.nd

shapiro.test(mbf$d21)       # n.nd
shapiro.test(mbfm$d21)      # nd
shapiro.test(mbff$d21)      # nd

library(gplots)
library(plotrix)

B6 <- "grey30"
BFMI <- "orange"
x <- seq(28,70,7)

#### Plot 1: Fat/Lean
# Fat Food, matB6
y <- c(mean(mb6$d28fat2lean), mean(mb6$d35fat2lean), mean(mb6$d42fat2lean), mean(mb6$d49fat2lean), mean(mb6$d56fat2lean), mean(mb6$d63fat2lean), mean(mb6$d70fat2lean))
ySD <- c(sd(mb6$d28fat2lean), sd(mb6$d35fat2lean), sd(mb6$d42fat2lean), sd(mb6$d49fat2lean), sd(mb6$d56fat2lean), sd(mb6$d63fat2lean), sd(mb6$d70fat2lean)) 

plot(x = c(28,70), y = c(0, 0.31), type='n', xlab="Age [days]", ylab="Fat Mass/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under HFD"))
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
points(x, y, cex=2, type='l', col=B6)
plotCI(x, y, uiw=ySD, add=T, gap=0, col=B6, pch=19, cex=2)

# Fat Food, matBFMI
y <- c(mean(mbf$d28fat2lean), mean(mbf$d35fat2lean), mean(mbf$d42fat2lean), mean(mbf$d49fat2lean), mean(mbf$d56fat2lean), mean(mbf$d63fat2lean), mean(mbf$d70fat2lean))
ySD <- c(sd(mbf$d28fat2lean), sd(mbf$d35fat2lean), sd(mbf$d42fat2lean), sd(mbf$d49fat2lean), sd(mbf$d56fat2lean), sd(mbf$d63fat2lean), sd(mbf$d70fat2lean)) 
plotCI(x, y, uiw=ySD, add=T, col=BFMI, gap=0, pch=16, cex=2)
points(x, y, cex=2, type='l', col=BFMI)
legend("bottomright", pch=c(16,16), col=c(B6,BFMI), c("m+f, matB6", "m+f, matBFMI"))
#### End of Plot 1: Fat/Lean

#### Plot 2: Bodyweight
#Fat Food, matB6 bodyweight
y <- c(mean(mb6$d28), mean(mb6$d35), mean(mb6$d42), mean(mb6$d49), mean(mb6$d56), mean(mb6$d63), mean(mb6$d70))
ySD <- c(sd(mb6$d28), sd(mb6$d35), sd(mb6$d42), sd(mb6$d49), sd(mb6$d56), sd(mb6$d63), sd(mb6$d70)) 

plot(x = c(28,70), y = c(0, 40), type='n', xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under HFD"))
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
points(x, y, cex=2, type='l', col=B6)
plotCI(x, y, uiw=ySD, add=T, gap=0, col=B6, pch=19, cex=2)

# Fat Food, matBFMI
y <- c(mean(mbf$d28), mean(mbf$d35), mean(mbf$d42), mean(mbf$d49), mean(mbf$d56), mean(mbf$d63), mean(mbf$d70))
ySD <- c(sd(mbf$d28), sd(mbf$d35), sd(mbf$d42), sd(mbf$d49), sd(mbf$d56), sd(mbf$d63), sd(mbf$d70)) 
plotCI(x, y, uiw=ySD, add=T, col=BFMI, gap=0, pch=16, cex=2)
points(x, y, cex=2, type='l', col=BFMI)
legend("bottomright", pch=c(16,16), col=c(B6,BFMI), c("m+f, matB6", "m+f, matBFMI"))
#### End of Plot 2: Bodyweight

#### Plot 3: Fat/Lean (males)
# FF, matB6, males
y <- c(mean(mb6m$d28fat2lean), mean(mb6m$d35fat2lean), mean(mb6m$d42fat2lean), mean(mb6m$d49fat2lean), mean(mb6m$d56fat2lean), mean(mb6m$d63fat2lean), mean(mb6m$d70fat2lean))
ySD <- c(sd(mb6m$d28fat2lean), sd(mb6m$d35fat2lean), sd(mb6m$d42fat2lean), sd(mb6m$d49fat2lean), sd(mb6m$d56fat2lean), sd(mb6m$d63fat2lean), sd(mb6m$d70fat2lean)) 

plot(x = c(28,70), y = c(0, 0.35), type='n', xlab="Age [days]", ylab="Fat/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under HFD"))
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
points(x, y, cex=2, type='l', col=B6)
plotCI(x, y, uiw=ySD, add=T, gap=0, col=B6, pch=19, cex=2)

# FF, matBFMI, males
y <- c(mean(mbfm$d28fat2lean), mean(mbfm$d35fat2lean), mean(mbfm$d42fat2lean), mean(mbfm$d49fat2lean), mean(mbfm$d56fat2lean), mean(mbfm$d63fat2lean), mean(mbfm$d70fat2lean))
ySD <- c(sd(mbfm$d28fat2lean), sd(mbfm$d35fat2lean), sd(mbfm$d42fat2lean), sd(mbfm$d49fat2lean), sd(mbfm$d56fat2lean), sd(mbfm$d63fat2lean), sd(mbfm$d70fat2lean)) 
plotCI(x, y, uiw=ySD, add=T, col=BFMI, gap=0, pch=16, cex=2)
points(x, y, cex=2, type='l', col=BFMI)
legend("bottomright", pch=c(15,15), col=c(B6,BFMI), c("matB6 males", "matBFMI males"))
#### End of Plot 3: Fat/Lean (males)

#### Plot 4: bodyweight (males)
# FF, matB6, males, 
y <- c(mean(mb6m$d28), mean(mb6m$d35), mean(mb6m$d42), mean(mb6m$d49), mean(mb6m$d56), mean(mb6m$d63), mean(mb6m$d70))
ySD <- c(sd(mb6m$d28), sd(mb6m$d35), sd(mb6m$d42n), sd(mb6m$d49), sd(mb6m$d56), sd(mb6m$d63), sd(mb6m$d70)) 

plot(x = c(28,70), y = c(0, 40), type='n', xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under HFD"))
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
points(x, y, cex=2, type='l', col=B6)
plotCI(x, y, uiw=ySD, add=T, gap=0, col=B6, pch=19, cex=2)

# FF, matBFMI, males
y <- c(mean(mbfm$d28), mean(mbfm$d35), mean(mbfm$d42), mean(mbfm$d49), mean(mbfm$d56), mean(mbfm$d63), mean(mbfm$d70))
ySD <- c(sd(mbfm$d28), sd(mbfm$d35), sd(mbfm$d42), sd(mbfm$d49), sd(mbfm$d56), sd(mbfm$d63), sd(mbfm$d70)) 
plotCI(x, y, uiw=ySD, add=T, col=BFMI, gap=0, pch=16, cex=2)
points(x, y, cex=2, type='l', col=BFMI)
legend("bottomright", pch=c(15,15), col=c(B6,BFMI), c("matB6 males", "matBFMI males"))
#### End of Plot 4: bodyweight (males)


#### Plot 5: Fat/Lean (females)
# FF, matB6, females
y <- c(mean(mb6f$d28fat2lean), mean(mb6f$d35fat2lean), mean(mb6f$d42fat2lean), mean(mb6f$d49fat2lean), mean(mb6f$d56fat2lean), mean(mb6f$d63fat2lean), mean(mb6f$d70fat2lean))
ySD <- c(sd(mb6f$d28fat2lean), sd(mb6f$d35fat2lean), sd(mb6f$d42fat2lean), sd(mb6f$d49fat2lean), sd(mb6f$d56fat2lean), sd(mb6f$d63fat2lean), sd(mb6f$d70fat2lean)) 

plot(x = c(28,70), y = c(0, 0.31), type='n', xlab="Age [days]", ylab="Fat/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under HFD"))
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
points(x, y, cex=2, type='l', col=B6)
plotCI(x, y, uiw=ySD, add=T, gap=0, col=B6, pch=19, cex=2)

# FF, matBFMI, females
y <- c(mean(mbff$d28fat2lean), mean(mbff$d35fat2lean), mean(mbff$d42fat2lean), mean(mbff$d49fat2lean), mean(mbff$d56fat2lean), mean(mbff$d63fat2lean), mean(mbff$d70fat2lean))
ySD <- c(sd(mbff$d28fat2lean), sd(mbff$d35fat2lean), sd(mbff$d42fat2lean), sd(mbff$d49fat2lean), sd(mbff$d56fat2lean), sd(mbff$d63fat2lean), sd(mbff$d70fat2lean)) 
plotCI(x, y, uiw=ySD, add=T, col=BFMI, gap=0, pch=16, cex=2)
points(x, y, cex=2, type='l', col=BFMI)
legend("bottomright", pch=c(16,16), col=c(B6,BFMI), c("females, matB6", "females, matBFMI"))
#### End of Plot 5: Fat/Lean (females)




#### stripchart fat/lean
plot(x = c(0,8), y = c(0, 0.4), cex.lab=1.5, ylab="Mass (g)", las=1, xlim=c(0,8), main=paste("F1 Population of Reciprocal B6N x BFMI860 Cross (HFD),\n", "Fat/Lean Ratio"), xlab="", xaxt='n', type="n")
grid(7, NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
stripchart(mb6$d28fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=0.2,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d35fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=1.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d42fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=2.55, col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d49fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=3.8,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d56fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=5.1,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d63fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=6.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d70fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=7.45, col=B6,    method="jitter", jitter = 0.05)
stripchart(mbf$d28fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=0.5,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d35fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=1.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d42fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=2.9,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d49fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=4.15, col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d56fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=5.4,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d63fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=6.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d70fat2lean, add=T, vertical=T, cex=1.2, pch=16, at=7.8,  col=BFMI,  method="jitter", jitter = 0.05)

text(c(0.3,1.5,2.8,4,5.2,6.5,7.7),-0.03 ,cex=1.2,srt=0, labels=seq(28,70,7), xpd=TRUE, col=1)
text(4,-0.06 ,cex=1.2,srt=0, labels="Age(days)", xpd=TRUE, col=1)
legend("topleft", inset=c(0.0,0.00), pch=c(16),col=c(B6,BFMI), c("maternal B6N","maternal BFMI860"), bg="white")

# segments draws a line
#FFB6   <- dataFFmatB6$totalGF
#FFBFMI <- dataFFmatBFMI$totalGF
#ZFB6   <- dataZFmatB6$totalGF
#ZFBFMI <- dataZFmatBFMI$totalGF
#segments(0.6, median(FFB6, na.rm=T), 1.0, (median(FFB6, na.rm=T)), col="grey", lwd=2)
#segments(0.6, median(FFBFMI, na.rm=T), 1.0, (median(FFBFMI, na.rm=T)), col="yellow3", lwd=2)
#segments(1.2, median(ZFB6, na.rm=T), 1.6, (median(ZFB6, na.rm=T)), col="grey", lwd=2)
#segments(1.2, median(ZFBFMI, na.rm=T), 1.6, (median(ZFBFMI, na.rm=T)), col="yellow3", lwd=2)

################# stripchart mrimass (FAT+LEAN)
# d28-d70 fat/lean
plot(x = c(0,8), y = c(10, 50), cex.lab=1.5, ylab="Mass (g)", las=1, xlim=c(0,8), main=paste("F1 Population of Reciprocal B6N x BFMI860 Cross (HFD),\n", "Fat+Lean"), xlab="", xaxt='n', type="n")
grid(7, NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
stripchart(mb6$d28mrimass, add=T, vertical=T, cex=1, pch=16, at=0.2,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d35mrimass, add=T, vertical=T, cex=1, pch=16, at=1.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d42mrimass, add=T, vertical=T, cex=1, pch=16, at=2.55, col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d49mrimass, add=T, vertical=T, cex=1, pch=16, at=3.8,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d56mrimass, add=T, vertical=T, cex=1, pch=16, at=5.1,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d63mrimass, add=T, vertical=T, cex=1, pch=16, at=6.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6$d70mrimass, add=T, vertical=T, cex=1, pch=16, at=7.45, col=B6,    method="jitter", jitter = 0.05)
stripchart(mbf$d28mrimass, add=T, vertical=T, cex=1, pch=16, at=0.5,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d35mrimass, add=T, vertical=T, cex=1, pch=16, at=1.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d42mrimass, add=T, vertical=T, cex=1, pch=16, at=2.9,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d49mrimass, add=T, vertical=T, cex=1, pch=16, at=4.15, col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d56mrimass, add=T, vertical=T, cex=1, pch=16, at=5.4,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d63mrimass, add=T, vertical=T, cex=1, pch=16, at=6.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbf$d70mrimass, add=T, vertical=T, cex=1, pch=16, at=7.8,  col=BFMI,  method="jitter", jitter = 0.05)

text(c(0.3,1.5,2.8,4,5.2,6.5,7.7),-3.5 ,cex=1,srt=0, labels=seq(28,70,7), xpd=TRUE, col=1)
text(4,-6 ,cex=1,srt=0, labels="Age(days)", xpd=TRUE, col=1)
legend("topleft", inset=c(0.0,0.00), pch=c(16),col=c(B6,BFMI), c("maternal B6N","maternal BFMI860"), bg="white")


################# stripchart mrimass (FAT+LEAN) m vs. f
# d28-d70 fat+lean, männchen
plot(x = c(0,8), y = c(10, 50), cex.lab=1.5, ylab="Mass (g)", las=1, xlim=c(0,8), main=paste("F1 Population of Reciprocal B6N x BFMI860 Cross (HFD),\n", "Fat+Lean"), xlab="", xaxt='n', type="n")
grid(7, NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
stripchart(mb6m$d28mrimass, add=T, vertical=T, cex=1, pch=15, at=0.2,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6m$d35mrimass, add=T, vertical=T, cex=1, pch=15, at=1.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6m$d42mrimass, add=T, vertical=T, cex=1, pch=15, at=2.55, col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6m$d49mrimass, add=T, vertical=T, cex=1, pch=15, at=3.8,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6m$d56mrimass, add=T, vertical=T, cex=1, pch=15, at=5.1,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6m$d63mrimass, add=T, vertical=T, cex=1, pch=15, at=6.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6m$d70mrimass, add=T, vertical=T, cex=1, pch=15, at=7.45, col=B6,    method="jitter", jitter = 0.05)
stripchart(mbfm$d28mrimass, add=T, vertical=T, cex=1, pch=15, at=0.5,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbfm$d35mrimass, add=T, vertical=T, cex=1, pch=15, at=1.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbfm$d42mrimass, add=T, vertical=T, cex=1, pch=15, at=2.9,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbfm$d49mrimass, add=T, vertical=T, cex=1, pch=15, at=4.15, col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbfm$d56mrimass, add=T, vertical=T, cex=1, pch=15, at=5.4,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbfm$d63mrimass, add=T, vertical=T, cex=1, pch=15, at=6.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbfm$d70mrimass, add=T, vertical=T, cex=1, pch=15, at=7.8,  col=BFMI,  method="jitter", jitter = 0.05)

# weibchen
stripchart(mb6f$d28mrimass, add=T, vertical=T, cex=1, pch=16, at=0.2,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6f$d35mrimass, add=T, vertical=T, cex=1, pch=16, at=1.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6f$d42mrimass, add=T, vertical=T, cex=1, pch=16, at=2.55, col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6f$d49mrimass, add=T, vertical=T, cex=1, pch=16, at=3.8,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6f$d56mrimass, add=T, vertical=T, cex=1, pch=16, at=5.1,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6f$d63mrimass, add=T, vertical=T, cex=1, pch=16, at=6.3,  col=B6,    method="jitter", jitter = 0.05)
stripchart(mb6f$d70mrimass, add=T, vertical=T, cex=1, pch=16, at=7.45, col=B6,    method="jitter", jitter = 0.05)
stripchart(mbff$d28mrimass, add=T, vertical=T, cex=1, pch=16, at=0.5,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbff$d35mrimass, add=T, vertical=T, cex=1, pch=16, at=1.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbff$d42mrimass, add=T, vertical=T, cex=1, pch=16, at=2.9,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbff$d49mrimass, add=T, vertical=T, cex=1, pch=16, at=4.15, col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbff$d56mrimass, add=T, vertical=T, cex=1, pch=16, at=5.4,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbff$d63mrimass, add=T, vertical=T, cex=1, pch=16, at=6.6,  col=BFMI,  method="jitter", jitter = 0.05)
stripchart(mbff$d70mrimass, add=T, vertical=T, cex=1, pch=16, at=7.8,  col=BFMI,  method="jitter", jitter = 0.05)

text(c(0.3,1.5,2.8,4,5.2,6.5,7.7),7 ,cex=1,srt=0, labels=seq(28,70,7), xpd=TRUE, col=1)
text(4,5 ,cex=1,srt=0, labels="Age(days)", xpd=TRUE, col=1)
legend("topleft", inset=c(0.0,0.00), pch=c(15,15,16,16),col=c(B6,BFMI), c("matB6N males","matBFMI860 males","matB6N females","matBFMI860 females"), bg="white")



####### statistics
f <- read.table(file="20140718_F1V.txt", sep="\t", na.strings = "NA", stringsAsFactors=F, header=TRUE )

anova(lm(f$totalGF ~ as.factor(f$WG) * f$sex + f$futter + f$dirocross)) 


anova(lm(f$d21 ~ as.factor(f$WG) * f$sex + f$futter + f$dirocross)) 
anova(lm(f$d28 ~ as.factor(f$WG) * f$sex + f$futter + f$dirocross))
model <- lm(f$d56fat2lean ~ as.factor(f$WG) + f$sex + f$futter)
model2 <- lm(model$residuals ~ f$dirocross)
anova(model2) # gruppieren WG in 2 gruppen 
boxplot(model$residuals ~ as.factor(f$dirocross))

anova(lm(f$glukosemg ~ as.factor(f$WG) * f$sex + f$futter + f$dirocross)) # gruppieren WG in 2 gruppen 

# m+f
anova(lm(f$d56fat2lean ~ as.factor(f$WG) * f$sex + f$futter + f$dirocross))
anova(lm(f$d63fat2lean ~ as.factor(f$WG) * f$sex + f$futter + f$dirocross))
anova(lm(f$d70fat2lean ~ as.factor(f$WG) * f$sex + f$futter + f$dirocross))
mean(f$d56fat2lean) #0.120
sd(f$d56fat2lean) #0.057


m <- subset(f , f$sex == "m")
f <- subset(f, f$sex=="f")
nrow(m)
nrow(f)

head(f)
mb6m <- subset(f, f$dirocross=="matB6" & sex == "m")
mbfm <- subset(f, f$dirocross=="matBFMI" & sex == "m")
nrow(mb6m)
nrow(mbfm)
mean(mb6m$d56fat2lean);sd(mb6m$d56fat2lean) # 0.0963 , 0.0368
mean(mbfm$d56fat2lean);sd(mbfm$d56fat2lean) # 0.1650 , 0.0813 # 13 tiere pro gruppe für power von 0.8

mb6 <- subset(f, f$dirocross=="matB6" )
mbf <- subset(f, f$dirocross=="matBFMI" )
nrow(mb6)
nrow(mbf)
mean(mb6$d56fat2lean);sd(mb6$d56fat2lean) # 0.1005 , 0.0331
mean(mbf$d56fat2lean);sd(mbf$d56fat2lean) # 0.1391 , 0.0698 # > getpower(powercurve, 0.80), You will need in each group: 30 animals



# m
anova(lm(m$d49fat2lean ~ as.factor(m$WG)  + m$futter + m$dirocross))
anova(lm(m$d56fat2lean ~ as.factor(m$WG)  + m$futter + m$dirocross))
anova(lm(m$d63fat2lean ~ as.factor(m$WG)  + m$futter + m$dirocross))
anova(lm(m$d70fat2lean ~ as.factor(m$WG)  + m$futter + m$dirocross))

mean(m$d56fat2lean)# 0.130
sd(m$d56fat2lean) #0.071

# f
anova(lm(f$d56fat2lean ~ as.factor(f$WG) + f$futter + f$dirocross))
anova(lm(f$d63fat2lean ~ as.factor(f$WG) + f$futter + f$dirocross))
anova(lm(f$d70fat2lean ~ as.factor(f$WG) + f$futter + f$dirocross))


infile <- "20140718_F1V.txt" 
f <- read.table(file=infile, sep="\t", na.strings = "NA", stringsAsFactors=F, header=TRUE )

fat <- c(f$d28mrifat, f$d35mrifat, f$d42mrifat, f$d49mrifat, f$d56mrifat, f$d63mrifat, f$d70mrifat)
lean <- c(f$d28mrilean, f$d35mrilean, f$d42mrilean, f$d49mrilean, f$d56mrilean, f$d63mrilean, f$d70mrilean)
fat2lean <- c(f$d28fat2lean, f$d35fat2lean, f$d42fat2lean, f$d49fat2lean, f$d56fat2lean, f$d63fat2lean, f$d70fat2lean)
bw <- c(f$d28, f$d35, f$d42, f$d49, f$d56, f$d63, f$d70)

age <- c(rep(28, 60), rep(35, 60), rep(42, 60), rep(49, 60), rep(56, 60), rep(63, 60), rep(70, 60))
wg <- rep(f$WG, 7)
sex <- rep(f$sex, 7)
futter<- rep(f$futter, 7)
dirocross<- rep(f$dirocross, 7)

anova(lm(bw ~ as.factor(age) + as.factor(wg) + sex + futter + dirocross))
modelFat <- lm(fat ~ as.factor(age) + as.factor(wg) + sex + futter + dirocross)
anova(modelFat)
modelLean <- lm(lean ~ as.factor(age) + as.factor(wg) + sex + futter + dirocross)
anova(modelLean)
anova(lm(fat2lean ~ as.factor(age) + as.factor(wg) + sex + futter + dirocross))


library(gplots)

boxplot(lean ~ as.factor(wg) + sex )
boxplot(lean ~ as.factor(wg))
boxplot(fat2lean ~ as.factor(wg))
boxplot(lean ~ sex + as.factor(wg))
anova(lm(lean ~ sex + as.factor(wg)))
anova(lm(lean ~ sex + as.factor(wg) + dirocross))
boxplot(lean ~ as.factor(wg) + dirocross+ sex, las=2)
b <- "grey50"
f <- "orange"
boxplot(lm(fat~age)$residuals ~ as.factor(wg) + dirocross+ sex, las=2, cex.axis=0.7, col = c(b,b,b,b,b,f,f,f,f,f), main="Residuals Fat Mass ~ age, for dirocross and litter size")
boxplot(lm(lean~age)$residuals ~ as.factor(wg) + dirocross+ sex, las=2, cex.axis=0.7,col = c(b,b,b,b,b,f,f,f,f,f), main="Residuals Lean Mass ~ age, for dirocross and litter size", varwidth=TRUE)




#oder so, boxplot2 macht das zwar automatisch, aber es ist unten halb verdeckt, also kleiner trick mit text, n bekommt man aber auch vom normalen boxplot
a <- boxplot2(lm(fat~age)$residuals ~ as.factor(wg) + dirocross+ sex, las=2, cex.axis=0.7, main="Residuals Fat Mass ~ age, for dirocross and litter size", shrink=0.6, top=F,textcolor="blue")
a <- a$n
datapoints <- paste( "i=", a, sep="")
datapoints 
boxplot(lm(fat~age)$residuals ~ as.factor(wg) + dirocross+ sex, las=2, cex.axis=0.7, col = c(b,b,b,b,b,f,f,f,f,f), main="Residuals Fat Mass ~ age, for dirocross and litter size")
#grid(NA, NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"),equilogs = FALSE)
abline(v = 0.5:20.5, lty = 1, lwd = 1, col = "grey50" )
abline(v = 10.5, lty = 1, lwd = 1, col = "grey10" )
text(1:20,-3.6, datapoints , cex=0.6, col="darkblue")
tiere <- lapply(a, "/", 7)
tiere <- unlist(tiere)
tiere <- paste("n=",tiere, sep="")
text(1:20,8.6, tiere, cex=0.6, col="red")
boxplot(lm(fat~age)$residuals ~ as.factor(wg) + dirocross+ sex, las=2, cex.axis=0.7, col = c(b,b,b,b,b,f,f,f,f,f), main="Residuals Fat Mass ~ age, for dirocross and litter size", add=T)



trait <- "Lean"
a <- boxplot2(lm(lean~age)$residuals ~ as.factor(wg) + dirocross+ sex, las=2, cex.axis=0.7, main="Residuals Fat Mass ~ age, for dirocross and litter size", shrink=0.6, top=F,textcolor="blue")
a <- a$n
n <- paste( "n=", a, sep="")
n
boxplot(lm(lean~age)$residuals ~ as.factor(wg) + dirocross+ sex, las=2, cex.axis=0.7, col = c(b,b,b,b,b,f,f,f,f,f), main=paste("Residuals ",trait," Mass ~ age, for dirocross and litter size"))
text(1:20,-6.7, n, cex=0.6, col="darkblue")

t.test(mb6m$d49,mbfm$d49)				# p-value = 
t.test(mb6m$d56,mbfm$d56)				# p-value = 
t.test(mb6m$d63,mbfm$d63)				# p-value = 
t.test(mb6m$d70,mbfm$d70)				# p-value = 

t.test(mb6f$d28,mbff$d28)
t.test(mb6f$d49,mbff$d49)				# p-value = 
t.test(mb6f$d56,mbff$d56)				# p-value = 
t.test(mb6f$d63,mbff$d63)				# p-value = 
t.test(mb6f$d70,mbff$d70)				# p-value = 


t.test(mb6m$d28mrifat,mbfm$d28mrifat)
t.test(mb6m$d49mrifat,mbfm$d49mrifat)		# p-value = 0.01784
t.test(mb6m$d56mrifat,mbfm$d56mrifat)		# p-value = 0.004927
t.test(mb6m$d63mrifat,mbfm$d63mrifat)		# p-value = 0.006829
t.test(mb6m$d70mrifat,mbfm$d70mrifat)		# p-value = 0.001552

t.test(mb6f$d28mrifat,mbff$d28mrifat)
t.test(mb6f$d49mrifat,mbff$d49mrifat)		# p-value = 
t.test(mb6f$d56mrifat,mbff$d56mrifat)		# p-value = 
t.test(mb6f$d63mrifat,mbff$d63mrifat)		# p-value = 
t.test(mb6f$d70mrifat,mbff$d70mrifat)		# p-value = 

t.test(mb6m$glukosemg,mbfm$glukosemg)		# p-value = 0.4022
t.test(mb6f$glukosemg,mbff$glukosemg)		# p-value = 0.5756

t.test(mb6m$totalGF,mbfm$totalGF)			# p-value = 0.001202
t.test(mb6f$totalGF,mbff$totalGF)			# p-value = 0.519

t.test(mb6m$totalRF,mbfm$totalRF)			# p-value = 0.0003214
t.test(mb6f$totalRF,mbff$totalRF)			# p-value = 0.7496

t.test(mb6m$muskel,mbfm$muskel)			# p-value = 0.6366
t.test(mb6f$muskel,mbff$muskel)			# p-value =  0.1005

t.test(mb6m$leber,mbfm$leber)				# p-value = 0.1278
t.test(mb6f$leber,mbff$leber)				# p-value = 0.3223

t.test(mb6m$BAT,mbfm$BAT)				# p-value = 0.0003848
t.test(mb6f$BAT,mbff$BAT)				# p-value = 0.9694

t.test(mb6m$LD,mbfm$LD)					# p-value = 0.9366
t.test(mb6f$LD,mbff$LD)					# p-value = 0.1845

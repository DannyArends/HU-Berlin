library(StAMPP)

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
# Read the genotypes formatted for STAMPP input
stammpinput <- t(read.table("StammpInput.txt", sep="\t"))
stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies

# Calculate pairwise FST values
stammpinput.fst <- stamppFst(stammpinput.freq, 10, 95, 4) # Population Fst values
stammpinput.fst$Fsts

#  Nei's genetic distance between individuals
stammp.D.pop <- stamppNeisD(stammpinput.freq, FALSE)

# Calculate AMOVA
stamppAmova(stammp.D.pop, stammpinput.freq, 10000)

library(emma)
ks <- c()  # holding the calculated kinships between 2 individuals
for(x in 1:100){
  #Generate 100 markers randomly for individual 1
  i1 <- round(runif(100))
  #Generate 100 markers randomly for individual 2
  i2 <- round(runif(100))
  ks <- c(ks,emma.kinship(cbind(i1,i2)) [1,2])
}
boxplot(ks)

ks <- c()  # holding the calculated kinships between 2 individuals
for(x in 1:100){
  #Generate 100 markers randomly for individual 1
  i1 <- round(runif(100))
  #change the first 50 markers randomly for individual 2
  i2 <- i1
  i2[1:50] <- round(runif(50))
  ks <- c(ks,emma.kinship(cbind(i1,i2)) [1,2])
}
boxplot(ks)


library(StAMPP)
p1 <- matrix(rep(rep("AA", 10), 2), 10, 2)
rownames(p1) <- paste0("P1_I", seq(1,10))
p1[,1] <- "AB"

p2 <- matrix(rep(rep("AB", 10), 2), 10, 2)
rownames(p2) <- paste0("P2_I", seq(1,10))
p2[,1] <- "BB"

p3 <- matrix(rep(rep("BB", 10), 2), 10, 2)
rownames(p3) <- paste0("P3_I", seq(1,10))
p3[,1] <- "AA"

combined <- rbind(p1,p2,p3)

exData <- cbind("Samples" = rownames(combined), Pop = c(rep("A", 10),rep("B", 10),rep("C", 10)),  Ploidy = 2, Format = "BiA", combined)
exData.freq <- stamppConvert(exData, "r") # Frequencies
exData.fst <- stamppFst(exData.freq, 1000) # Population Fst values
exData.fst$Fsts

stammp.D.pop <- stamppNeisD(exData.freq, FALSE)    #  Nei's genetic distance between individuals
stamppAmova(stammp.D.pop, exData.freq, 10000)      # Calculate AMOVA
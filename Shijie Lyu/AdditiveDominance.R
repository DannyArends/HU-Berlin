# Example to show / learn about additive genetic effects and dominance deviations

#genotypes <- as.character(unlist(read.table("D:/Sjije/AddDomData.txt", row.names=1)))
#phenotype <- as.numeric(unlist(read.table("D:/Sjije/AddDomData2.txt", header=TRUE,colClasses="character")))

nind <- 1000
nperm <- 50

effRange <- seq(0.01, 0.5, 0.01)

matri <- NULL

for(percAA in seq(0.01, 0.49, 0.01)){
  # AA, AB, BB genotypes frequencies, we set hetrozygous to be 50 % of the samples
  percentage <- c(percAA , 0.50); percentage <- c(percentage, 1-sum(percentage))

  sumofsq <- NULL
  for(effect in effRange){
    # Simulation parameters
    sim.dom <- -effect
    sim.add <- 1
    genotypes <- c(rep("AA", ceiling(percentage[1] * nind)), rep("AB", ceiling(percentage[2] * nind)), rep("BB", ceiling(percentage[3] * nind)))
    means <- NULL
    coefs <- NULL
    for(p in 1:nperm){
      marker <- genotypes[sample(nind,nind)]
      markerAdditive  <- as.numeric(as.factor(marker))
      markerDominance <- as.numeric(marker == "AB")

      phenotype <- rnorm(nind, 10, 0.5)
      
      # Additive effect A, H, B 
      #phenotype[markerAdditive==1] <- phenotype[markerAdditive==1] + 0 * sim.add
      #phenotype[markerAdditive==2] <- phenotype[markerAdditive==2] + 1 * sim.add
      #phenotype[markerAdditive==3] <- phenotype[markerAdditive==3] + 2 * sim.add

      # Major dominance / hybrid vigor (heterosis) effect of genotype H
      phenotype[markerAdditive==2] <- phenotype[markerAdditive==2] + sim.dom

      means <- rbind(means, c(Additive = (mean(phenotype[markerAdditive==3]) - mean(phenotype[markerAdditive==1]))/2, Dominance = (mean(phenotype[markerAdditive==2]) - mean(phenotype[markerAdditive!=2])) ))
      m012  <- markerAdditive
      coefs <- rbind(coefs, lm(phenotype ~  markerDominance + markerAdditive)$coefficients)
    }
    sumofsq <- c(sumofsq, sum((means[,"Dominance"] - coefs[,"markerDominance"])^2))
    cat("Done", percAA, " ", effect,"\n")
  }
  matri <- rbind(matri, (sumofsq/effRange))
}

colz <- colorRampPalette(c(rgb(1,0,0,1), rgb(0,1,0,1), rgb(0,0,1,1)), alpha = TRUE)(nrow(matri))

plot(x=c(0,dim(matri)[2]), y = c(0,max(matri)))
for(i in 1:nrow(matri)){
  points(matri[i,], col = colz[i], t='l')
}

image(t(log(1 + matri)), col=colorRampPalette(c(rgb(1,1,1,1), rgb(1,1,0,1), rgb(1,0,0,1)), alpha = TRUE)(8),breaks = c(0, 0.01, 0.1, 0.5, 1, 2, 3, 10, 100),xaxt='n',yaxt='n')
axis(1, at=seq(0.0, 0.9, 0.01), effRange)
axis(2, at=seq(0.015, 0.975, 0.04), seq(0.01, 0.49, 0.02),las=2)
grid(length(effRange), length(seq(0.01, 0.49, 0.01)), col="white",lty=1)

#plot(c(0,2), c(0,2), xlab = "(Mean(aa) - Mean(AA)) / 2", ylab="Linear regression coefficient", t = 'n')
#points(means[,"Additive"], coefs[,"markerAdditive"], pch=18)
#cor(means, coefs, method="spearman")
#plot(means[,"Dominance"], coefs[,"markerDominance"], pch=18)
#plot(sumofsq / seq(0, 10, 0.2),t='l')
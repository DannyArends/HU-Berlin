# Example to show / learn about additive genetic effects and dominance deviations

nind <- 1000
nperm <- 10000

# AA, AB, BB genotypes frequencies
percentage <- c(0.25 , 0.50); percentage <- c(percentage, 1-sum(percentage))

# Simulation parameters
sim.dom <- -30
sim.add <- 1

genotypes <- c(rep("AA", percentage[1] * nind), rep("AB", percentage[2] * nind), rep("BB", percentage[3] * nind))

means <- NULL
coefs <- NULL
for(x in 1:nperm){
  marker <- genotypes[sample(nind,nind)]
  markerAdditive  <- as.numeric(as.factor(marker))
  markerDominance <- as.numeric(marker == "AB")

  phenotype <- rnorm(nind, 10, 0.5)
  
  # Additive effect A, H, B 
  phenotype[markerAdditive==1] <- phenotype[markerAdditive==1] + 0 * sim.add
  phenotype[markerAdditive==2] <- phenotype[markerAdditive==2] + 1 * sim.add
  phenotype[markerAdditive==3] <- phenotype[markerAdditive==3] + 2 * sim.add

  # Major dominance / hybrid viggor effect of genotype H
  phenotype[markerAdditive==2] <- phenotype[markerAdditive==2] + sim.dom

  means <- rbind(means, c(Additive = (mean(phenotype[markerAdditive==3]) - mean(phenotype[markerAdditive==1]))/2, Dominance = (mean(phenotype[markerAdditive==2]) - mean(phenotype[markerAdditive!=2])) ))
  m012  <- markerAdditive
  coefs <- rbind(coefs, lm(phenotype ~  markerDominance + markerAdditive)$coefficients)
}

plot(c(0,2), c(0,2), xlab = "(Mean(aa) - Mean(AA)) / 2", ylab="Linear regression coefficient", t = 'n')
points(means[,"Additive"], coefs[,"markerAdditive"], pch=18)
cor(means, coefs, method="spearman")

plot(means[,"Dominance"], coefs[,"markerDominance"], pch=18)
#
# Simulation for DFG project Gudrun
#

# Generate a genetic marker at HWE with a given minor allele frequency
generateMarker <- function(MAF = 0.1) {
  marker <- rep(1, 100)
  P <- MAF
  Q = 1 - MAF
  nA <- round((P ^ 2) * 100,0)
  nH <- round(2 * (P * Q) * 100,0)
  nB <- round((Q ^ 2) * 100,0)
  marker[sample(1:100, nH)] <-  0
  while(nA > 0){
    idx <- sample(1:100, 1)
    if (marker[idx] == 1) {
      marker[idx] <- -1
      nA <- nA - 1
    }
  }
  return(marker)
}

# Given a genetic marker, generate a phenotype where the marker has a given effect size
generatePhenotype <- function(marker, effectsize = 0.4){
  phenotype <- rnorm(100, 8, 2)
  m <- mean(phenotype)
  phenotype[marker == 1] <- phenotype[marker == 1] + effectsize * m
  phenotype[marker == 0] <- phenotype[marker == 0] + (0.5 * m * effectsize)
  return(phenotype)
}

ntests <- 100 # Number of test to perform to compute the power
niter <- 20 # Number of times the power calculation is repeated
nmarkers <- c(25, 50, 100, 500, 1000) # Number of markers in upstream region
mafs <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5) # Minor allele frequencies of markers tested
effectsizes <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0)

res <- c()
for(nmar in nmarkers){
  for(maf in mafs){
    mym <- generateMarker(maf)
    for(e in effectsizes){
      for(i in 1:niter){
        nD <- 0
        for(x in 1:ntests){
          ph <- generatePhenotype(mym, e)
          pval <- anova(lm(ph ~ mym))[[5]][1]
          if( pval < (0.05 / nmar) ){
            nD = nD + 1
          }
        }
        res <- rbind(res, c(i, nmar, maf, e, nD / ntests))
      }
    }
    cat("Done maf", maf, "\n")
  }
  cat("Done marker", nmar, "\n")
}
colnames(res) <- c("iteration", "nmar", "maf", "eff", "power")

# Plot the results of the simulation
op <- par(mfrow = c(3, 2))
colorz <- colorRampPalette(c("red", "green"))(length(effectsizes))
for(nmar in nmarkers){
  toplot <- res[res[,"nmar"] == nmar,]
  plot(x = c(0, 0.5), y = c(0, 100), t = 'n', xlab = 'Minor allele frequency', ylab = "Power", las = 2, sub = paste0("# of markers in 1Mb upstream: ", nmar))
  n <- 1
  for(x in unique(toplot[,"eff"])){
    medians <- c()
    for(m in unique(toplot[,"maf"])){
      toplott <- toplot[which(toplot[, "eff"] == x & toplot[, "maf"] == m),]
      boxplot(100 * toplott[, "power"], at = m, add = TRUE, boxwex = 0.015, col = colorz[n], yaxt='n')
      medians <- c(medians, median(100 * toplott[, "power"]))
    }
    points(x = c(0, unique(toplot[,"maf"])), c(0,medians), t = 'l', col = colorz[n])
    n <- n + 1
  }
  abline(h = 80, col = "green")
}
plot(c(0,1), c(0,1), t = 'n', xaxt='n', yaxt='n', xlab="", ylab="", bty='n')
legend("toplef", paste0("Effect size: ", round(effectsizes,1)), fill=colorz)

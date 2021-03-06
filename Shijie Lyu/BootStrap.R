# Bootstrapping a QTL confidence interval
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2015
# first written Jun, 2015

gS <- round(runif(200))
genotypes <- NULL
for(x in 1:100){
  genotypes <- cbind(genotypes, gS)
  gS[round(runif(1,1,100))] <- gS[round(runif(1,1,100))]
  gS[round(runif(1,1,100))] <- gS[round(runif(1,1,100))]
  gS[round(runif(1,1,100))] <- gS[round(runif(1,1,100))]
  gS[round(runif(1,1,100))] <- gS[round(runif(1,1,100))]
}

phenotypes <- rnorm(200) + 0.5 * genotypes[,15]

mapQTL <- function(genotypes,phenotypes){
  pvalues <- NULL
  for(x in 1:ncol(genotypes)){
    pvalues <- c(pvalues, anova(lm(phenotypes ~ genotypes[,x]))[[5]][1])
  }
  pvalues
}

QTLs <- -log10(mapQTL(genotypes, phenotypes))

boots <- NULL
for(x in 1:100){
  idx <- sample(nrow(genotypes), 200, replace=TRUE)
  boots <- rbind(boots, mapQTL(genotypes[idx, ], phenotypes[idx]))
}
top <- apply(-log10(boots), 2, function(x){ sort(x)[97] })
bot <- apply(-log10(boots), 2, function(x){ sort(x)[2] })

perms <- NULL
for(x in 1:100){
  perms <- rbind(perms, mapQTL(genotypes[sample(nrow(genotypes)), ], phenotypes[idx]))
}
threshold <- apply(-log10(perms), 2, function(x){ sort(x)[97] })
maxes <- apply(-log10(perms), 1, function(x){ max(x) })

plot(c(0,100), c(0,45), t='n')
points(QTLs, t='l', lwd=2)
points(top, t = 'l', col='red')
points(top, t = 'h', col='red')
points(QTLs ,t='h', col="blue")
points(bot, t = 'l', col='blue')
points(bot ,t='h', col="white")
points(threshold, t='l', lwd=2, lty=3)
abline(h=sort(maxes)[95], lwd=2,lty=3)

abline(h=bot[which.max(QTLs)], lwd=2,lty=2)

which(QTLs > bot[which.max(QTLs)])
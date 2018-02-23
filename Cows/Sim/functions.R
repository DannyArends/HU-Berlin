# Calculate the effects of the different markers
doMarkers <- function(genotypes, markers, effect.type = c("add", "dom")){
  change <- rep(0, ncol(genotypes))
  for(marker in names(markers)){
    he <- which(genotypes[marker,] == "0/1")
    ho <- which(genotypes[marker,] == "1/1")
    change[he] = change[he] + markers[marker]
    if(effect.type[1] == "add") change[ho] = change[ho] + 2 * markers[marker]
    if(effect.type[1] == "dom") change[ho] = change[ho] + markers[marker]
  }
  return(change)
}

# Calculate alternative allele frequency
calc.allele.freq <- function(gts){
  gt <- table(unlist(strsplit(as.character(gts), "/")))
  freq <- (gt / sum(gt))
  return(freq["1"])
}

simulateLactation <- function(genotypes, production.markers, disease.markers, production.per.lactation = 10000, natural.variance = 0.10, doPlot=TRUE){
  cows.before <- ncol(genotypes)
  disease.base <- round(rnorm(ncol(genotypes), 100, 1), 1)
  disease.mod <- doMarkers(genotypes, disease.markers)
  disease.chance <- disease.base + disease.mod
  D <- round(runif(length(disease.chance), min(disease.chance), max(disease.chance) + 1), 1)

  sick <- which(D > disease.chance)
  sick <- sample(sick, (round(rnorm(1, 30, 1), 1)/100) * cows.before)   # ~ 30 % of the herd has some health issues
  health <- sample(sick, (round(rnorm(1, 20, 1), 1)/100) * cows.before) # ~ 30 % of the herd is removed because of health reasons

  # Uwe: on average ~ 7.500 to 10.000 Liters per cow per lacatation
  production.base <- round(rnorm(ncol(genotypes), production.per.lactation, natural.variance * production.per.lactation), 1)
  production.mod <- doMarkers(genotypes, production.markers)
  production.lactation <- production.base + production.mod

  production.lactation[sick] = production.lactation[sick] * .75
  
  production.sorted.ind <- sort(production.lactation, index.return=T)$ix
  production.sorted.ind <- production.sorted.ind[-which(production.sorted.ind %in% health)]
  production.sorted.ind <- production.sorted.ind[-which(production.sorted.ind %in% sick)]
  culled <- production.sorted.ind[1:round((round(rnorm(1, 10, 0.5), 1)/100) * length(production.sorted.ind))]

  genotypes <- genotypes[, -unique(c(culled,health))]
  if (doPlot) {
    mcol <- rep("black", cows.before)
    mcol[sick] <- "orange"
    mcol[culled] <- "red"
    mcol[health] <- "blue"
    par(mfrow=c(1, 2))
    plot(disease.chance, D, pch=19, col = mcol, xlab = "Disease resistance", ylab="Disease risk encountered (D)")
    plot(1:length(production.lactation), production.lactation, xlab="Individual", ylab="Simulated milk yield", pch=19, col = mcol)
  }

  #Uwe: ~30 % of animals are removed from the population per lactation
  #Uwe: ~5 to 10 % of animals are removed because of leistungs problems, The rest due to health issues
  perc.sick <- paste0("(", round(length(sick)/cows.before * 100, 1),"%)")
  perc.health <- paste0("(", round(length(health)/cows.before * 100, 1),"%)")
  perc.culled <- paste0("(", round(length(culled)/cows.before * 100, 1),"%)")
  cat(lactation, "Before:", cows.before, "Sick:", length(sick), perc.sick, "Health:", length(health),perc.health, "Culled:", length(culled),perc.culled, "Left:", ncol(genotypes), "\n")

  return(genotypes)
}

plotFreq.change <- function(freqs.change, map, chromosomes, chr.starts, chr.lengths, production.markers, disease.markers) {
  mcol = rep("black", nrow(map))
  mcol[rownames(map) %in% names(disease.markers)] <- "firebrick1"
  mcol[rownames(map) %in% names(production.markers)] <- "black"

  mpch = rep("", nrow(map))
  mpch[rownames(map) %in% names(which(disease.markers < 0))] <- "-"
  mpch[rownames(map) %in% names(which(disease.markers > 0))] <- "+"
  mpch[rownames(map) %in% names(which(production.markers > 0))] <- "+"
  mpch[rownames(map) %in% names(which(production.markers < 0))] <- "-"

  chrID <- 1 + as.numeric(factor(map[,"CHROM"], levels=chromosomes)) %% 2

  par(mfrow=c(1, 1))
  plot(c(0, max(as.numeric(map[,"CPOS"]))), c(-0.1, 0.1), xlab="Chromosome", ylab="Allele frequency change", xaxt='n', t='n', las=2)
  points(as.numeric(map[,"CPOS"]), freqs.change, col=c("gray80", "#C6DBEF")[chrID], pch=19, cex=0.5, xlab="Chromosome", xaxt='n')
  points(as.numeric(map[,"CPOS"]), freqs.change, col=mcol, pch=mpch, cex=2)
  axis(1, at = chr.starts[1:length(chromosomes)] + (chr.lengths)/2, gsub("Chr", "", chromosomes))
}


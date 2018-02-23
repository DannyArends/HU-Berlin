setwd("D:/Edrive/Cow/Fugato")

genotypes <- read.csv("FugateNsequencing_genotypes.txt", sep = "\t", colClasses = "character", check.names = FALSE)
phedata <- read.csv("FugateNsequencing_samples.txt", sep = "\t", colClasses = "character", check.names = FALSE)
annotation <- read.csv("FugateNsequencing_probes.txt", sep = "\t", colClasses = "character", check.names = FALSE)

setwd("D:/Ddrive/Papers/idea Uwe/Sim")
source("functions.R")

chromosomes <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10")

HF <- which(phedata[,"breed"] == "HF")
map <- annotation[,c("CHROM", "POS")]
map <- map[which(map[,"CHROM"] %in% chromosomes),]
map <- map[order( map[,"CHROM"], as.numeric(map[,"POS"])), ]
map <- map[which(rownames(map) %in% rownames(genotypes)), ]

genotypes <- genotypes[rownames(map), HF]

call.rates <- apply(genotypes,1,function(x){1 - (sum(is.na(x)) / length(x))})
freqs.before <- apply(genotypes, 1, calc.allele.freq)

maf <- freqs.before
maf[which(maf > 0.5)] <- 1 - maf[which(maf > 0.5)]

allele.freq.low <- which(freqs.before <= 0.05 | freqs.before >= 0.95)
call.rate.low <- which(call.rates < 0.95)

to.remove <- unique(c(call.rate.low,allele.freq.low)) #Markers to remove because of missing values / low allele frequency

map <- map[-to.remove,]
genotypes <- genotypes[-to.remove,]
freqs.before <- freqs.before[-to.remove]

map <- cbind(map, "CPOS" = 0)
chr.starts <- c(0)
chr.lengths <- c()
lengthsofar <- 0
gap <- 10000000
for(chr in chromosomes){
  onChr <- which(map[,"CHROM"] == chr)
  chrL <- max(as.numeric(map[onChr,"POS"]))
  map[onChr, "CPOS"] <- as.numeric(map[onChr, "POS"]) + lengthsofar
  lengthsofar <- lengthsofar + chrL + gap
  chr.starts <- c(chr.starts, lengthsofar)
  chr.lengths <- c(chr.lengths, chrL)
}

#Variables, we might want to investigate
set.seed(0)
source("functions.R")
production.per.lactation <- 9000
natural.variance <- 0.05          # 10 % of the milk production is random

n.simulations <- 100              # Number of simulations to perform
n.lactations <- 4                 # Number of lactations simulated per simulation

n.production.markers <- 5         # Number of genetic markers affecting production
n.disease.markers <- 10           # Number of genetic markers affecting disease

production.data <- vector("list", n.simulations)
disease.data <- vector("list", n.simulations)
for(simulation in 1:n.simulations) {
  cat("---- Simulation:", simulation, "----\n")
  genodata <- genotypes

  # Uwe: Production markers show a change ~ 0.1 and 1.5 % of leistung
  production.markers <- round(runif(n.production.markers, 0, 0.05 * production.per.lactation), 0) * sign(runif(n.production.markers) - 0.5)
  names(production.markers) <- sample(rownames(map), n.production.markers)
  production.markers <- (production.markers - mean(production.markers))

  # Uwe: Perhaps separate out into difference diseases (Health issues, Leistungs affecting)
  disease.markers <- round(runif(n.disease.markers, 1, 5), 1) * sign(runif(n.disease.markers) - 0.5)
  names(disease.markers) <- sample(rownames(map), n.disease.markers)
  disease.markers <- (disease.markers - mean(disease.markers))

  for(lactation in 1:n.lactations) {
    genodata <- simulateLactation(genodata, production.markers, disease.markers, production.per.lactation, natural.variance)
  }

  freqs.after <- apply(genodata, 1, calc.allele.freq)
  freqs.change <- (freqs.after - freqs.before)

  plotFreq.change(freqs.change, map, chromosomes, chr.starts, chr.lengths, production.markers, disease.markers)

  production.stats <- cbind(production.change = production.markers, freq.before = freqs.before[names(production.markers)], freq.change = freqs.change[names(production.markers)])
  disease.stats <- cbind(disease.chance = disease.markers, freq.before = freqs.before[names(disease.markers)], freq.change = freqs.change[names(disease.markers)])
  
  production.data[[simulation]] <- list(mean.change = mean(freqs.change), sd.change=sd(freqs.change), production.stats = production.stats)
  disease.data[[simulation]] <- list(mean.change = mean(freqs.change), sd.change=sd(freqs.change), disease.stats = disease.stats)
}

production.effectsize <- (unlist(lapply(production.data, function(x){x$production.stats[,1]})) / production.per.lactation) * 100
production.sdchange <- unlist(lapply(production.data, function(x){x$production.stats[,3] / x$sd.change}))
production.model <- lm(production.sdchange ~ production.effectsize+0)
production.plottext <- paste0("Regression ß: ", round(production.model$coefficients["production.effectsize"],2))

disease.effectsize <- unlist(lapply(disease.data, function(x){x$disease.stats[,1]}))
disease.sdchange <- unlist(lapply(disease.data, function(x){x$disease.stats[,3] / x$sd.change}))
disease.model <- lm(disease.sdchange ~ disease.effectsize+0)
disease.plottext <- paste0("Regression ß: ", round(disease.model$coefficients["disease.effectsize"],2))

par(mfrow=c(1, 2))
plot(production.effectsize, production.sdchange, main="Production", sub = production.plottext, xlab = "effectsize (%)", ylab="# of SDs AlleleFreq Change", pch=19, mgp=c(2,1,0))
abline(h = c(2,-2), col = "orange", lty=2)
abline(h = c(3,-3), col = "chartreuse3", lty=2)
abline(v = c(1.5, -1.5), lty=2)
abline(a = 0, b = production.model$coefficients["production.effectsize"], col="coral3")
legend("topleft", legend=c("Regression line", "± 2 SD", "± 3 SD", "1.5 % effectsize"), lwd=c(1, 1, 1, 1), col=c("coral3", "orange", "chartreuse3", "black"), lty=c(1, 2, 2, 2), bg="white")

plot(disease.effectsize, disease.sdchange, main="Disease", sub = disease.plottext, xlab = "effectsize (%)", ylab="# of SDs AlleleFreq Change", pch=19, mgp=c(2,1,0))
abline(h = c(2,-2), col = "orange", lty=2)
abline(h = c(3,-3), col = "chartreuse3", lty=2)
abline(a = 0, b = disease.model$coefficients["disease.effectsize"], col="coral3")
legend("topleft", legend=c("Regression line", "± 2 SD", "± 3 SD"), lwd=c(1, 1, 1), col=c("coral3", "orange", "chartreuse3"), lty=c(1, 2, 2), bg="white")

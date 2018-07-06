library(lme4)
setwd("D:/Edrive/Mouse/DNA/MegaMuga/inputF2")
genotypes <- read.table(file="cleaned_recoded_genotypes_F2.txt", sep = '\t', check.names = FALSE)
map <- read.table(file="cleaned_map_25MbGap.txt", sep = '\t')
phenotypes <- read.table(file="cleaned_phenotypes_F2.txt", sep = '\t')

# Extract covariates
timepoints <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")
daysSinceStart <- as.numeric(gsub("d", "", timepoints))-21
# Group everything into 1 big model
phenotype <- NULL
individual <- NULL
subfamily <- NULL
littersize <- NULL
litternumber <- NULL
season <- NULL
timepoint <- NULL
for(tp in timepoints){
  phenotype <- c(phenotype, phenotypes[, tp])
  individual <- c(individual, rownames(phenotypes))
  subfamily <- c(subfamily, phenotypes[, "Vater"])
  littersize <- c(littersize, phenotypes[, "WG2"])
  litternumber <- c(litternumber, as.character(phenotypes[, "W.Label"]))
  season <- c(season, as.character(phenotypes[, "Season"]))
  timepoint <- c(timepoint, rep(tp, nrow(phenotypes)))
}
timepoint <- as.numeric(gsub("d", "", timepoint)) - 21
litternumber5 <- as.factor(litternumber) # Fixed effect: Number of litter (factor)
litternumber2 <- as.factor(as.numeric(litternumber5 == "A") + 1)
subfamily <- as.factor(subfamily) # Fixed effect: Number of litter (factor)
littersize <- as.factor(littersize) # Fixed effect: Number of litter (factor)
topmarker <- factor(unlist(genotypes["UNC5048297", individual]), levels = c("N", "H", "B"))
names(subfamily) <- individual
ltype5 <- factor(paste0(litternumber5, "_", littersize), levels = c("A_8", "A_9", "A_10", "A_11", "A_12", "B_8", "B_9", "B_12", "C_8", "C_9", "D_8", "D_9","E_8", "E_9"))
ltype2 <- factor(paste0(litternumber2, "_", littersize), levels = c("1_8", "1_9", "1_12", "2_8", "2_9", "2_10", "2_11", "2_12"))
names(ltype5) <- individual
names(ltype2) <- individual
season <- as.factor(season)
names(season) <- individual

ctrl = lmerControl(optimizer ="Nelder_Mead", optCtrl=list(maxfun=20000))



plot.data <- function(...){
  set.seed(1)
  plot(c(0,50), c(0,60), t = 'n', xlab="Time since start (Day)", ylab="Bodyweight (Gram)", ...)
  for(day in daysSinceStart){
    onTP <- which(timepoint == day)
    points(x = rep(day, length(onTP)) + 0.4 * runif(length(onTP)), phenotype[onTP], col="black", pch=19)
  }
}

m0 <- lmer(phenotype ~ (1|individual), REML = FALSE, control = ctrl)
plot.data(main="P = μ + (1|I)")
eff.fixed.intercept <- fixef(m0)["(Intercept)"]
for(ind in rownames(phenotypes)){
  eff.random.intercept <- ranef(m0)$individual[ind, "(Intercept)"]
  
  curve(eff.fixed.intercept + 0 * x, add=TRUE, col=rgb(0.2,0.2,0.2,0.2))
}

m1 <- lmer(phenotype ~ timepoint + (1|individual), REML = FALSE, control = ctrl)
plot.data(main="P = μ + T + (1|I)")
eff.fixed.intercept <- fixef(m1)["(Intercept)"]
for(ind in rownames(phenotypes)){
  eff.fixed.timepoint <- fixef(m1)["timepoint"]
  
  eff.random.intercept <- ranef(m1)$individual[ind, "(Intercept)"]
  
  curve(eff.fixed.intercept + eff.fixed.timepoint * x, add=TRUE, col=rgb(0.2,0.2,0.2,0.2))
}

m2 <- lmer(phenotype ~ timepoint + I(timepoint^2) + (1|individual), REML = FALSE, control = ctrl)
plot.data(main="P = μ + T + T2 + (1|I)")
eff.fixed.intercept <- fixef(m2)["(Intercept)"]
for(ind in rownames(phenotypes)){
  eff.fixed.timepoint <- fixef(m2)["timepoint"]
  eff.fixed.timepoint2 <- fixef(m2)["I(timepoint^2)"]
  
  eff.random.intercept <- ranef(m2)$individual[ind, "(Intercept)"]
  
  curve(eff.fixed.intercept + eff.fixed.timepoint * x + eff.fixed.timepoint2 * x^2, add=TRUE, col=rgb(0.2,0.2,0.2,0.2))
}


m3 <- lmer(phenotype ~ timepoint + I(timepoint^2) + I(timepoint^3)  + (1|individual), REML = FALSE, control = ctrl)
plot.data(main="P = μ + T + T2 + T3 + (1|I)")
eff.fixed.intercept <- fixef(m3)["(Intercept)"]
for(ind in rownames(phenotypes)){
  eff.fixed.timepoint <- fixef(m3)["timepoint"]
  eff.fixed.timepoint2 <- fixef(m3)["I(timepoint^2)"]
  eff.fixed.timepoint3 <- fixef(m3)["I(timepoint^3)"]

  eff.random.intercept <- ranef(m3)$individual[ind, "(Intercept)"]

  curve(eff.fixed.intercept + eff.fixed.timepoint * x + eff.fixed.timepoint2 * x^2 + eff.fixed.timepoint3 * x^3, add=TRUE, col=rgb(0.2,0.2,0.2,0.2))
}

m4 <- lmer(phenotype ~ subfamily + timepoint + I(timepoint^2) + I(timepoint^3) + (1|individual), REML = FALSE, control = ctrl)
plot.data(main="P = μ + F + T + T2 + T3 + (1|I)")
eff.fixed.intercept <- fixef(m4)["(Intercept)"]
for(ind in rownames(phenotypes)){
  ind.subfamily <- as.character(subfamily[ind])
  eff.subfamily <- fixef(m4)[paste0("subfamily", ind.subfamily)]
  if(is.na(eff.subfamily)) eff.subfamily <- 0
  
  eff.fixed.timepoint <- fixef(m4)["timepoint"]
  eff.fixed.timepoint2 <- fixef(m4)["I(timepoint^2)"]
  eff.fixed.timepoint3 <- fixef(m4)["I(timepoint^3)"]

  eff.random.intercept <- ranef(m4)$individual[ind, "(Intercept)"]

  curve(eff.fixed.intercept +  eff.subfamily + + eff.fixed.timepoint * x + eff.fixed.timepoint2 * x^2 + eff.fixed.timepoint3 * x^3, add=TRUE, col=rgb(0.2,0.2,0.2,0.2))
}


m5 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + (1|individual), REML = FALSE, control = ctrl)
plot.data(main="P = μ + F + LT2 + T + T2 + T3 + (1|I)")
eff.fixed.intercept <- fixef(m5)["(Intercept)"]
for(ind in rownames(phenotypes)){
  ind.subfamily <- as.character(subfamily[ind])
  eff.subfamily <- fixef(m5)[paste0("subfamily", ind.subfamily)]
  if(is.na(eff.subfamily)) eff.subfamily <- 0
  
  ind.ltype <- as.character(ltype2[ind])
  eff.ltype <- fixef(m5)[paste0("ltype", ind.ltype)]
  if(is.na(eff.ltype)) eff.ltype <- 0

  eff.fixed.timepoint <- fixef(m5)["timepoint"]
  eff.fixed.timepoint2 <- fixef(m5)["I(timepoint^2)"]
  eff.fixed.timepoint3 <- fixef(m5)["I(timepoint^3)"]

  eff.random.intercept <- ranef(m5)$individual[ind, "(Intercept)"]

  curve(eff.fixed.intercept +  eff.subfamily + eff.ltype + eff.fixed.timepoint * x + eff.fixed.timepoint2 * x^2 + eff.fixed.timepoint3 * x^3, add=TRUE, col=rgb(0.2,0.2,0.2,0.2))
}


m6 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + (1|individual), REML = FALSE, control = ctrl)
plot.data(main="P = μ + F + LT2 + T + T2 + T3 + jObes1 + jObes1:T + (1|I)")
eff.fixed.intercept <- fixef(m6)["(Intercept)"]
for(ind in rownames(phenotypes)){
  ind.subfamily <- as.character(subfamily[ind])
  eff.subfamily <- fixef(m6)[paste0("subfamily", ind.subfamily)]
  if(is.na(eff.subfamily)) eff.subfamily <- 0
  
  ind.ltype <- as.character(ltype2[ind])
  eff.ltype <- fixef(m6)[paste0("ltype", ind.ltype)]
  if(is.na(eff.ltype)) eff.ltype <- 0

  eff.fixed.timepoint <- fixef(m6)["timepoint"]
  eff.fixed.timepoint2 <- fixef(m6)["I(timepoint^2)"]
  eff.fixed.timepoint3 <- fixef(m6)["I(timepoint^3)"]
  
  ind.topmarker <- as.character(topmarker[ind])
  eff.topmarker <- fixef(m6)[paste0("topmarker", ind.topmarker)]
  if (is.na(eff.topmarker)) eff.topmarker <- 0
  if (is.na(ind.topmarker)) eff.topmarker <- NA

  eff.topmarkertime <- fixef(m6)[paste0("timepoint:topmarker", ind.topmarker)]
  if (is.na(eff.topmarkertime)) eff.topmarkertime <- 0

  eff.random.intercept <- ranef(m6)$individual[ind, "(Intercept)"]

  if(ind.topmarker == "N") colz <- "blue"
  if(ind.topmarker == "H") colz <- "gray"
  if(ind.topmarker == "B") colz <- "orange"
  
  curve(eff.fixed.intercept +  eff.subfamily + eff.ltype + eff.topmarker + eff.topmarkertime * x + eff.fixed.timepoint * x + eff.fixed.timepoint2 * x^2 + eff.fixed.timepoint3 * x^3, add=TRUE, col=colz)
}
legend("topleft", c("jObes1 = B6N", "jObes1 = H", "jObes1 = BFMI"), col= c("blue","gray", "orange"), lwd=1)


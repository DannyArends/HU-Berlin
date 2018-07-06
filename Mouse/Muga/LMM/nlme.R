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


m0     <- lmer(phenotype ~ subfamily + (1|individual), REML = FALSE, control = ctrl)
m1_L2  <- lmer(phenotype ~ subfamily + litternumber2 + (1|individual), REML = FALSE, control = ctrl)
m1_L5  <- lmer(phenotype ~ subfamily + litternumber5 + (1|individual), REML = FALSE, control = ctrl)

m2_L2  <- lmer(phenotype ~ subfamily + litternumber2 + littersize + (1|individual), REML = FALSE, control = ctrl)
m2_L5  <- lmer(phenotype ~ subfamily + litternumber5 + littersize + (1|individual), REML = FALSE, control = ctrl)

m2_Lt2 <- lmer(phenotype ~ subfamily + ltype2 + (1|individual), REML = FALSE, control = ctrl)
m2_Lt5 <- lmer(phenotype ~ subfamily + ltype5 + (1|individual), REML = FALSE, control = ctrl)

# Data for the AIC table
diff(AIC(m0, m1_L2)[,"AIC"]) 
diff(AIC(m0, m1_L5)[,"AIC"])
diff(AIC(m0, m2_L2)[,"AIC"])
diff(AIC(m0, m2_L5)[,"AIC"])
diff(AIC(m0, m2_Lt2)[,"AIC"])
diff(AIC(m0, m2_Lt5)[,"AIC"])

diff(AIC(m1_L2, m1_L5)[,"AIC"])
diff(AIC(m1_L2, m2_L2)[,"AIC"])
diff(AIC(m1_L2, m2_L5)[,"AIC"])
diff(AIC(m1_L2, m2_Lt2)[,"AIC"])
diff(AIC(m1_L2, m2_Lt5)[,"AIC"])

diff(AIC(m1_L5, m2_L2)[,"AIC"])
diff(AIC(m1_L5, m2_L2)[,"AIC"])
diff(AIC(m1_L5, m2_Lt2)[,"AIC"])
diff(AIC(m1_L5, m2_Lt5)[,"AIC"])

diff(AIC(m2_L2, m2_L5)[,"AIC"])
diff(AIC(m2_L2, m2_Lt2)[,"AIC"])
diff(AIC(m2_L2, m2_Lt5)[,"AIC"])

diff(AIC(m2_L5, m2_Lt2)[,"AIC"])
diff(AIC(m2_L5, m2_Lt5)[,"AIC"])

diff(AIC(m2_Lt2, m2_Lt5)[,"AIC"])

# Data for the big LMM models
m0 <- lmer(phenotype ~ subfamily + ltype2 + (1|individual), REML = FALSE, control = ctrl)
m1 <- lmer(phenotype ~ subfamily + ltype2 + season + (1|individual), REML = FALSE, control = ctrl)
diff(AIC(m0, m1)[,"AIC"])
m2 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + (1|individual), REML = FALSE, control = ctrl)
diff(AIC(m1, m2)[,"AIC"])
m3 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
diff(AIC(m2, m3)[,"AIC"])
m4 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + (timepoint|individual), REML = FALSE, control = ctrl)
diff(AIC(m3, m4)[,"AIC"])
m5 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual), REML = FALSE, control = ctrl)
diff(AIC(m4, m5)[,"AIC"])
m6 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + I(timepoint^4) + (timepoint|individual), REML = FALSE, control = ctrl)
diff(AIC(m5, m6)[,"AIC"])
m7 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
diff(AIC(m5, m7)[,"AIC"])

top.markers <- c("UNC1938399","UNC030576333","JAX00522656","UNC9857021","UNC10016752","UNC090485124","UNC30294194")

library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)["R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  timepoints <- grep("timepoint", names(variance.explained))
  interactions <- grep(":", names(variance.explained))
  if(length(interactions) > 0) timepoints <- timepoints[-which(timepoints %in% interactions)]
  scaled <- variance.explained[-timepoints] / sum(variance.explained[-timepoints])
  return(scaled)
}



null <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual), REML = FALSE, control = ctrl)
variance.explained.raw(null)
variance.explained.scaled(null)

mnR0 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + (timepoint|individual), REML = TRUE, control = ctrl)
mnR0.raw <- variance.explained.raw(mnR0)
100 * (variance.explained.scaled(null)["Unexplained"] - variance.explained.scaled(mnR0)["Unexplained"])

nR1top <- factor(unlist(genotypes["UNC1938399", individual]), levels = c("N", "H", "B"))
mnR1 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + nR1top + nR1top:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
mnR1.raw <- variance.explained.raw(mnR1)
100 * (variance.explained.scaled(mnR0)["Unexplained"] - variance.explained.scaled(mnR1)["Unexplained"])


nR2top <- factor(unlist(genotypes["UNC030576333", individual]), levels = c("N", "H", "B"))
mnR2 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + nR2top + nR2top:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
mnR2.raw <- variance.explained.raw(mnR2)
100 * (variance.explained.scaled(mnR0)["Unexplained"] - variance.explained.scaled(mnR2)["Unexplained"])


nR3top <- factor(unlist(genotypes["JAX00522656", individual]), levels = c("N", "H", "B"))
mnR3 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + nR3top + nR3top:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
mnR3.raw <- variance.explained.raw(mnR3)
100 * (variance.explained.scaled(mnR0)["Unexplained"] - variance.explained.scaled(mnR3)["Unexplained"])


nR4top <- factor(unlist(genotypes["UNC9857021", individual]), levels = c("N", "H", "B"))
mnR4 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + nR4top + nR4top:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
mnR4.raw <- variance.explained.raw(mnR4)
100 * (variance.explained.scaled(mnR0)["Unexplained"] - variance.explained.scaled(mnR4)["Unexplained"])


nR5top <- factor(unlist(genotypes["UNC10016752", individual]), levels = c("N", "H", "B"))
mnR5 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + nR5top + nR5top:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
mnR5.raw <- variance.explained.raw(mnR5)
100 * (variance.explained.scaled(mnR0)["Unexplained"] - variance.explained.scaled(mnR5)["Unexplained"])


nR6top <- factor(unlist(genotypes["UNC090485124", individual]), levels = c("N", "H", "B"))
mnR6 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + nR6top + nR6top:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
mnR6.raw <- variance.explained.raw(mnR6)
100 * (variance.explained.scaled(mnR0)["Unexplained"] - variance.explained.scaled(mnR6)["Unexplained"])


nR7top <- factor(unlist(genotypes["UNC30294194", individual]), levels = c("N", "H", "B"))
mnR7 <- lmer(phenotype ~ subfamily + ltype2 + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + nR7top + nR7top:timepoint + (timepoint|individual), REML = FALSE, control = ctrl)
mnR7.raw <- variance.explained.raw(mnR7)
100 * (variance.explained.scaled(mnR0)["Unexplained"] - variance.explained.scaled(mnR7)["Unexplained"])

createEffectMatrix <- function(lmemodel, topmarker = NULL, marker = NULL, markername = NULL){
  fixedeffects <- fixef(lmemodel)
  randomeffects <- ranef(lmemodel)$individual
  effectmatrix <- NULL
  for(ind in rownames(phenotypes)){
    eff.fixed.intercept <- fixedeffects["(Intercept)"]

    ind.subfamily <- as.character(subfamily[ind])
    eff.subfamily <- fixedeffects[paste0("subfamily", ind.subfamily)]
    if (is.na(eff.subfamily)) eff.subfamily <- 0

    ind.ltype2 <- as.character(ltype2[ind])
    eff.ltype2 <- fixedeffects[paste0("ltype2", ind.ltype2)]
    if (is.na(eff.ltype2)) eff.ltype2 <- 0

    eff.fixed.timepoint <- fixedeffects["timepoint"]
    eff.fixed.timepoint2 <- fixedeffects["I(timepoint^2)"]
    eff.fixed.timepoint3 <- fixedeffects["I(timepoint^3)"]

    if (!is.null(topmarker)) {  # Jobes1 marker effects
      ind.topmarker <- as.character(topmarker[ind])
      eff.topmarker <- fixedeffects[paste0("topmarker", ind.topmarker)]
      if (is.na(eff.topmarker)) eff.topmarker <- 0
      if (is.na(ind.topmarker)) eff.topmarker <- NA

      eff.topmarkertime <- fixedeffects[paste0("timepoint:topmarker", ind.topmarker)]
      if (is.na(eff.topmarkertime)) eff.topmarkertime <- 0
    }
    if (!is.null(marker)) {
      ind.marker <- as.character(marker[ind])
      eff.marker <- fixedeffects[paste0(markername, ind.marker)]
      if (is.na(eff.marker)) eff.marker <- 0
      if (is.na(ind.marker)) eff.marker <- NA

      eff.markertime <- fixedeffects[paste0("timepoint:", markername, ind.marker)]
      if (is.na(eff.markertime)) eff.markertime <- 0
    }
    eff.random.intercept <- randomeffects[ind, "(Intercept)"]
    eff.random.timepoint <- randomeffects[ind, "timepoint"]
    
    all.effects <- c("intercept" = as.numeric(eff.fixed.intercept),
                     "subfamily" = as.numeric(eff.subfamily),
                     "ltype2" = as.numeric(eff.ltype2),
                     "timepoint" = as.numeric(eff.fixed.timepoint),
                     "I(timepoint^2)" = as.numeric(eff.fixed.timepoint2),
                     "I(timepoint^3)" = as.numeric(eff.fixed.timepoint3))
    if (!is.null(topmarker)) {  # Jobes1 marker effects
      all.effects <- c(all.effects, "topmarker" = as.numeric(eff.topmarker),
                                    "timepoint:topmarker" = as.numeric(eff.topmarkertime))
    }
    if (!is.null(marker)) {
      all.effects <- c(all.effects, "marker" = as.numeric(eff.marker),
                                    "timepoint:marker" = as.numeric(eff.markertime))
    }
    all.effects <- c(all.effects, "random:intercept" = as.numeric(eff.random.intercept),
                                  "random:timepoint" = as.numeric(eff.random.timepoint))
                                          
    effectmatrix <- rbind(effectmatrix, all.effects)
  }
  rownames(effectmatrix) <- rownames(phenotypes)
  
  # make sure that we remove any individuals which have NA estimates, since they will cause issues when predicting the expected value
  containsNA <- which(apply(effectmatrix,1,function(x){ any(is.na(x))}))
  if(length(containsNA) > 0) effectmatrix <- effectmatrix[-containsNA, ]
  return(effectmatrix)
}

getCorSqEstimates <- function(observed, predicted) {
  cors <- NULL
  for(day in as.character(daysSinceStart)){
    cors <- c(cors, cor(predicted[,day], observed[rownames(predicted), day], method="pearson") ^ 2)
  }
  return(cors)
}


# Observed phenotype values at daysSinceStart
observed <- phenotypes[,timepoints]
colnames(observed) <- daysSinceStart

# Calculate effects and predicted values for the null model
effects.null <- createEffectMatrix(null)
effects.null[1:5, ]
estimates.null <- apply(effects.null, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3, from=0, to=50, add=FALSE)
})
predicted.null <- matrix(unlist(lapply(estimates.null, "[", 'y')), nrow(effects.null), length(estimates.null[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.null), estimates.null[[1]]$x))
predicted.null <- predicted.null[,colnames(observed)]
predicted.null[1:5, ]
varExp.null <- getCorSqEstimates(observed, predicted.null)

# Calculate effects and predicted values for the mnR0 model
effects.mnR0 <- createEffectMatrix(mnR0, topmarker)
effects.mnR0[1:5, ]
estimates.mnR0 <- apply(effects.mnR0, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR0 <- matrix(unlist(lapply(estimates.mnR0, "[", 'y')), nrow(effects.mnR0), length(estimates.mnR0[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR0), estimates.mnR0[[1]]$x))
predicted.mnR0 <- predicted.mnR0[,colnames(observed)]
predicted.mnR0[1:5, ]
varExp.mnR0 <- getCorSqEstimates(observed, predicted.mnR0)

# Calculate effects and predicted values for the mnR1 model
effects.mnR1 <- createEffectMatrix(mnR1, topmarker, nR1top, "nR1top")
effects.mnR1[1:5, ]
estimates.mnR1 <- apply(effects.mnR1, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x + eff["marker"] + eff["timepoint:marker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR1 <- matrix(unlist(lapply(estimates.mnR1, "[", 'y')), nrow(effects.mnR1), length(estimates.mnR1[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR1), estimates.mnR1[[1]]$x))
predicted.mnR1 <- predicted.mnR1[,colnames(observed)]
predicted.mnR1[1:5, ]
varExp.mnR1 <- getCorSqEstimates(observed, predicted.mnR1)

# Calculate effects and predicted values for the mnR2 model
effects.mnR2 <- createEffectMatrix(mnR2, topmarker, nR2top, "nR2top")
effects.mnR2[1:5, ]
estimates.mnR2 <- apply(effects.mnR2, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x + eff["marker"] + eff["timepoint:marker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR2 <- matrix(unlist(lapply(estimates.mnR2, "[", 'y')), nrow(effects.mnR2), length(estimates.mnR2[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR2), estimates.mnR2[[1]]$x))
predicted.mnR2 <- predicted.mnR2[,colnames(observed)]
predicted.mnR2[1:5, ]
varExp.mnR2 <- getCorSqEstimates(observed, predicted.mnR2)

# Calculate effects and predicted values for the mnR3 model
effects.mnR3 <- createEffectMatrix(mnR3, topmarker, nR3top, "nR3top")
effects.mnR3[1:5, ]
estimates.mnR3 <- apply(effects.mnR3, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x + eff["marker"] + eff["timepoint:marker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR3 <- matrix(unlist(lapply(estimates.mnR3, "[", 'y')), nrow(effects.mnR3), length(estimates.mnR3[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR3), estimates.mnR3[[1]]$x))
predicted.mnR3 <- predicted.mnR3[,colnames(observed)]
predicted.mnR3[1:5, ]
varExp.mnR3 <- getCorSqEstimates(observed, predicted.mnR3)

# Calculate effects and predicted values for the mnR4 model
effects.mnR4 <- createEffectMatrix(mnR4, topmarker, nR4top, "nR4top")
effects.mnR4[1:5, ]
estimates.mnR4 <- apply(effects.mnR4, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x + eff["marker"] + eff["timepoint:marker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR4 <- matrix(unlist(lapply(estimates.mnR4, "[", 'y')), nrow(effects.mnR4), length(estimates.mnR4[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR4), estimates.mnR4[[1]]$x))
predicted.mnR4 <- predicted.mnR4[,colnames(observed)]
predicted.mnR4[1:5, ]
varExp.mnR4 <- getCorSqEstimates(observed, predicted.mnR4)

# Calculate effects and predicted values for the mnR5 model
effects.mnR5 <- createEffectMatrix(mnR5, topmarker, nR5top, "nR5top")
effects.mnR5[1:5, ]
estimates.mnR5 <- apply(effects.mnR5, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x + eff["marker"] + eff["timepoint:marker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR5 <- matrix(unlist(lapply(estimates.mnR5, "[", 'y')), nrow(effects.mnR5), length(estimates.mnR5[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR5), estimates.mnR5[[1]]$x))
predicted.mnR5 <- predicted.mnR5[,colnames(observed)]
predicted.mnR5[1:5, ]
varExp.mnR5 <- getCorSqEstimates(observed, predicted.mnR5)

# Calculate effects and predicted values for the mnR6 model
effects.mnR6 <- createEffectMatrix(mnR6, topmarker, nR6top, "nR6top")
effects.mnR6[1:5, ]
estimates.mnR6 <- apply(effects.mnR6, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x + eff["marker"] + eff["timepoint:marker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR6 <- matrix(unlist(lapply(estimates.mnR6, "[", 'y')), nrow(effects.mnR6), length(estimates.mnR6[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR6), estimates.mnR6[[1]]$x))
predicted.mnR6 <- predicted.mnR6[,colnames(observed)]
predicted.mnR6[1:5, ]
varExp.mnR6 <- getCorSqEstimates(observed, predicted.mnR6)

# Calculate effects and predicted values for the mnR7 model
effects.mnR7 <- createEffectMatrix(mnR7, topmarker, nR7top, "nR7top")
effects.mnR7[1:5, ]
estimates.mnR7 <- apply(effects.mnR7, 1, function(eff){
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + eff["timepoint"]*x + eff["I(timepoint^2)"]*x^2 + eff["I(timepoint^3)"]*x^3 + 
        eff["topmarker"] + eff["timepoint:topmarker"]*x + eff["marker"] + eff["timepoint:marker"]*x, from=0, to=50, add=FALSE)
})
predicted.mnR7 <- matrix(unlist(lapply(estimates.mnR7, "[", 'y')), nrow(effects.mnR7), length(estimates.mnR7[[1]]$x), byrow=TRUE, dimnames=list(rownames(effects.mnR7), estimates.mnR7[[1]]$x))
predicted.mnR7 <- predicted.mnR7[,colnames(observed)]
predicted.mnR7[1:5, ]
varExp.mnR7 <- getCorSqEstimates(observed, predicted.mnR7)

varExplained <- cbind(jObes1 = varExp.mnR0 - varExp.null, 
                      mnR1 = varExp.mnR1 - varExp.mnR0, 
                      mnR2 = varExp.mnR2 - varExp.mnR0,
                      mnR3 = varExp.mnR3 - varExp.mnR0,
                      mnR4 = varExp.mnR4 - varExp.mnR0,
                      mnR5 = varExp.mnR5 - varExp.mnR0,
                      mnR6 = varExp.mnR6 - varExp.mnR0,
                      mnR7 = varExp.mnR7 - varExp.mnR0) * 100
rownames(varExplained) <- daysSinceStart
varExplained <- round(varExplained, 2)

write.table(varExplained, "varExplained.txt",sep = "\t")



### OLD code, there be dragons here



op <- par(mfrow=c(1,1))
mcol <- 2
atAdj = -0.7
plot(c(0,50), c(0,60), t = 'n', main = "LMM estimates, no random effect per individual", xlab="Time since start (Day)", ylab="Bodyweight (Gram)")
for(type in as.character(unique(topmarker))) {
  for(day in daysSinceStart){
    onTP <- which(timepoint == day & topmarker == type)
    points(x = rep(day,length(onTP)) + atAdj + 0.2 * runif(length(onTP)), phenotype[onTP], col=mcol)
  }
  atAdj = atAdj + 0.7
  mcol = mcol + 1
}

colorz <- c(rgb(0.1,0.1,0.8,0.2), rgb(0.1,0.8,0.1,0.1), rgb(0.8,0.1,0.1,0.1))

cnt <- 0
estimates <- apply(effectmatrix, 1, function(eff){
  cnt <<- cnt + 1
  ind <- rownames(effectmatrix)[cnt]
  curve(eff["intercept"] + eff["subfamily"] + eff["ltype2"] + 
        eff["timepoint"] * x + eff["I(timepoint^2)"] * x^2 + eff["I(timepoint^3)"]  * x^3,
        #eff["topmarker"] + eff["timepoint:topmarker"] * x, #+ eff["random:intercept"] + eff["random:timepoint"] * x, 
        add=TRUE, col=colorz[topmarker[ind]])
})
observed <- phenotypes[,timepoints]
colnames(observed) <- daysSinceStart

predicted <- matrix(unlist(lapply(estimates, "[", 'y')), nrow(effectmatrix), length(estimates[[1]]$x), byrow=TRUE, dimnames=list(rownames(effectmatrix), estimates[[1]]$x))
predicted <- predicted[,colnames(observed)]

residual <- observed - predicted

cors <- NULL
for(day in as.character(daysSinceStart)){
  cors <- c(cors, cor(predicted[,day], observed[,day], method="spearman") ^ 2)
  cat(day, " ", cor(predicted[,day], observed[,day], method="spearman") ^ 2, "\n")
}



sumOfSquares <- function(x){ return(sum(((x - mean(x))^2))) }




for(ind in rownames(phenotypes)){
  ind.subfamily <- as.character(subfamily[ind])
  eff.subfamily <- lmemodel$coefficients$fixed[paste0("subfamily", ind.subfamily)]
  if(is.na(eff.subfamily)) eff.subfamily <- 0

  ind.ltype <- as.character(ltype[ind])
  eff.ltype <- lmemodel$coefficients$fixed[paste0("ltype", ind.ltype)]
  if(is.na(eff.ltype)) eff.ltype <- 0

  ind.season <- as.character(season[ind])
  eff.season <- lmemodel$coefficients$fixed[paste0("season", ind.season)]
  if(is.na(eff.season)) eff.season <- 0

  ind.topmarker <- as.character(topmarker[ind])
  eff.topmarker <- lmemodel$coefficients$fixed[paste0("topmarker", ind.topmarker)]
  if(is.na(eff.topmarker)) eff.topmarker <- 0

  eff.markertime <- lmemodel$coefficients$fixed[paste0("timepoint:topmarker", ind.topmarker)]
  if(is.na(eff.markertime)) eff.markertime <- 0

  eff.rnd.intercept <- lmemodel$coefficients$random$individual[ind, "(Intercept)"]
  eff.rnd.timepoint <- lmemodel$coefficients$random$individual[ind, "timepoint"]
  
  curve(lmemodel$coefficients$fixed["(Intercept)"] +
        eff.subfamily +
        eff.ltype +
        eff.season +
        lmemodel$coefficients$fixed["timepoint"] * x +
        lmemodel$coefficients$fixed["I(timepoint^2)"] * x^2 + 
        lmemodel$coefficients$fixed["I(timepoint^3)"] * x^3 +
        eff.topmarker +
        eff.markertime * x
        , add=TRUE, col=c(rgb(0.8,0.1,0.1,0.1), rgb(0.1,0.8,0.1,0.1), rgb(0.1,0.1,0.8,0.2))[topmarker[ind]])

  #cat(ind, " ", ind.subfamily, "=", eff.subfamily, " ", ind.ltype, "=", eff.ltype, " ", ind.season, "=", eff.season, "\n")
}


mcol <- 2
atAdj = -0.7
plot(c(0,70), c(0, 60), t = 'n', main = "LMM estimates, random effects per individual", xlab="Time (Day)", ylab="Bodyweight (Gram)")
for(type in as.character(unique(topmarker))) {
  for(tp in as.numeric(gsub("d", "", timepoints))){
    onTP <- which(timepoint == tp & topmarker == type)
    points(x = rep(tp,length(onTP)) + atAdj + 0.2 * runif(length(onTP)), phenotype[onTP], col=mcol)
  }
  atAdj = atAdj + 0.7
  mcol = mcol + 1
}

for(ind in rownames(phenotypes)){
  ind.subfamily <- as.character(subfamily[ind])
  eff.subfamily <- lmemodel$coefficients$fixed[paste0("subfamily", ind.subfamily)]
  if(is.na(eff.subfamily)) eff.subfamily <- 0

  ind.ltype <- as.character(ltype[ind])
  eff.ltype <- lmemodel$coefficients$fixed[paste0("ltype", ind.ltype)]
  if(is.na(eff.ltype)) eff.ltype <- 0

  ind.season <- as.character(season[ind])
  eff.season <- lmemodel$coefficients$fixed[paste0("season", ind.season)]
  if(is.na(eff.season)) eff.season <- 0

  ind.topmarker <- as.character(topmarker[ind])
  eff.topmarker <- lmemodel$coefficients$fixed[paste0("topmarker", ind.topmarker)]
  if(is.na(eff.topmarker)) eff.topmarker <- 0

  eff.markertime <- lmemodel$coefficients$fixed[paste0("timepoint:topmarker", ind.topmarker)]
  if(is.na(eff.markertime)) eff.markertime <- 0

  eff.rnd.intercept <- lmemodel$coefficients$random$individual[ind, "(Intercept)"]
  eff.rnd.timepoint <- lmemodel$coefficients$random$individual[ind, "timepoint"]
  
  curve(lmemodel$coefficients$fixed["(Intercept)"] + eff.rnd.intercept +
        eff.subfamily +
        eff.ltype +
        eff.season +
        lmemodel$coefficients$fixed["timepoint"] * x + eff.rnd.timepoint * x +
        lmemodel$coefficients$fixed["I(timepoint^2)"] * x^2 + 
        lmemodel$coefficients$fixed["I(timepoint^3)"] * x^3 +
        eff.topmarker +
        eff.markertime * x
        , add=TRUE, col=c(rgb(0.8,0.1,0.1,0.1), rgb(0.1,0.8,0.1,0.1), rgb(0.1,0.1,0.8,0.2))[topmarker[ind]])

  #cat(ind, " ", ind.subfamily, "=", eff.subfamily, " ", ind.ltype, "=", eff.ltype, " ", ind.season, "=", eff.season, "\n")
}


# Chromosome variables
chrs <- 1:20
names(chrs) <- c(1:19, "X")

chrs.starts <- c(0)
chrs.lengths <- c()
chrs.summed <- 0
chr.gap <- 25000000

for(chr in names(chrs)){
  onChr <- which(map[,"Chr"] == chr)
  chr.length <- max(as.numeric(map[onChr, "Mb_NCBI38"]))
  chrs.summed <- chrs.summed + chr.length + chr.gap
  chrs.lengths <- c(chrs.lengths, chr.length + chr.gap)
  chrs.starts <- c(chrs.starts, chrs.summed)
}
chrs.lengths <- c(chrs.lengths, NA)


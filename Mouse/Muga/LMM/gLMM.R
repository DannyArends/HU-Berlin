# LMM MQM mapping of time series data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
library(parallel)

setwd("~/LMMMQM")
genotypes <- read.table(file="cleaned_recoded_genotypes_F2.txt", sep = '\t', check.names = FALSE)
map <- read.table(file="cleaned_map_25MbGap.txt", sep = '\t')
phenotypes <- read.table(file="cleaned_phenotypes_F2.txt", sep = '\t')

# Extract covariates
timepoints <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")
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
litternumber <- as.factor(as.numeric(litternumber == "A") + 1)
ltype <- factor(paste0(litternumber, "_", littersize), levels = c("1_8", "1_9", "1_12", "2_8", "2_9", "2_10", "2_11", "2_12"))
names(ltype) <- individual
subfamily <- as.factor(subfamily)
topmarker <- factor(unlist(genotypes["UNC5048297", individual]), levels = c("N", "H", "B"))
season <- as.factor(season)
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

### Map everything using an across timepoint linear mixed model, using subfamily as a random effect and individual as grouping factor
cl <- makeCluster(getOption("cl.cores", 48))
clusterExport(cl, "phenotype") # export phenotype to the nodes
clusterExport(cl, "individual") # export individual to the nodes
clusterExport(cl, "subfamily") # export subfamily to the nodes
clusterExport(cl, "ltype") # export ltype to the nodes
clusterExport(cl, "topmarker") # export topmarker to the nodes
clusterExport(cl, "timepoint") # export timepoint to the nodes
clusterEvalQ(cl, library(lme4)) # load lme4 per node
LMMglobal <- parApply(cl, genotypes[, individual], 1, function(marker){
  ctrl = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun=20000))
  noG <- which(is.na(marker))
  MF <- factor(as.character(marker), levels = c("N", "H", "B"))
  if(length(noG) > 0){
    phenotype <- phenotype[-noG]
    individual <- individual[-noG]
    subfamily <- subfamily[-noG]
    ltype <- ltype[-noG]
    timepoint <- timepoint[-noG]
    MF <- MF[-noG]
    topmarker <- topmarker[-noG]
  }

  pvals <- rep(NA, 12)
  names(pvals) <- c("F_M_C", "M_N_C", "F_N_C", "F_M", "M_N", "F_N", "timepoint:MFB", "timepoint:MFH", "timepoint:MFN", "timepoint:MFCB", "timepoint:MFCH", "timepoint:MFCN")
  tryCatch(
    lmemodel.fullC <- lmer(phenotype ~ subfamily + ltype + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + MF + MF:timepoint + (timepoint|individual), control=ctrl, REML = FALSE)
    , error = function(e) e)
  tryCatch(
    lmemodel.markerC <- lmer(phenotype ~ subfamily + ltype + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + MF + (timepoint|individual), control=ctrl, REML = FALSE)
    , error = function(e) e)
  tryCatch(
    lmemodel.nullC <- lmer(phenotype ~ subfamily + ltype + timepoint + I(timepoint^2) + I(timepoint^3) + topmarker + topmarker:timepoint + (timepoint|individual), control=ctrl, REML = FALSE)
    , error = function(e) e)
  tryCatch(
    lmemodel.full <- lmer(phenotype ~ subfamily + ltype + timepoint + I(timepoint^2) + I(timepoint^3) + MF + MF:timepoint + (timepoint|individual), control=ctrl, REML = FALSE)
    , error = function(e) e)
  tryCatch(
    lmemodel.marker <- lmer(phenotype ~ subfamily + ltype + timepoint + I(timepoint^2) + I(timepoint^3) + MF + (timepoint|individual), control=ctrl, REML = FALSE)
    , error = function(e) e)
  tryCatch(
    lmemodel.null <- lmer(phenotype ~ subfamily + ltype + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual), control=ctrl, REML = FALSE)
    , error = function(e) e)
    
  tryCatch(pvals["F_M_C"] <- anova(lmemodel.fullC, lmemodel.markerC)["lmemodel.fullC", "Pr(>Chisq)"], error = function(e) e)
  tryCatch(pvals["M_N_C"] <- anova(lmemodel.markerC, lmemodel.nullC)["lmemodel.markerC", "Pr(>Chisq)"], error = function(e) e)
  tryCatch(pvals["F_N_C"] <- anova(lmemodel.fullC, lmemodel.nullC)["lmemodel.fullC", "Pr(>Chisq)"], error = function(e) e)
  tryCatch(pvals["F_M"] <- anova(lmemodel.full, lmemodel.marker)["lmemodel.full", "Pr(>Chisq)"], error = function(e) e)
  tryCatch(pvals["M_N"] <- anova(lmemodel.marker, lmemodel.null)["lmemodel.marker", "Pr(>Chisq)"], error = function(e) e)
  tryCatch(pvals["F_N"] <- anova(lmemodel.full, lmemodel.null)["lmemodel.full", "Pr(>Chisq)"], error = function(e) e)
  tryCatch(pvals[7:9] <- summary(lmemodel.full)$coefficients[, "Estimate"][c("timepoint:MFB", "timepoint:MFH", "timepoint:MFN")], error = function(e) e)
  tryCatch(pvals[10:12] <- summary(lmemodel.fullC)$coefficients[, "Estimate"][c("timepoint:MFB", "timepoint:MFH", "timepoint:MFN")], error = function(e) e)
  pvals[names(which(is.na(pvals[7:12])))] <- 0
  mG <- names(which(table(MF) == 0))
  if(length(mG) > 0){
    pvals[paste0("timepoint:MF",mG)] <- NA
    pvals[paste0("timepoint:MFC",mG)] <- NA
  }
  return(pvals)
})
stopCluster(cl)
rownames(LMMglobal) <- c("F_M_C", "M_N_C", "F_N_C", "F_M", "M_N", "F_N", "timepoint:MFB", "timepoint:MFH", "timepoint:MFN", "timepoint:MFCB", "timepoint:MFCH", "timepoint:MFCN")

write.table(t(LMMglobal), file=paste0("LMM_Global.txt"), sep='\t', quote=FALSE)


results.lm <- vector("list", length(timepoints))
names(results.lm) <- timepoints
results.lmm <- vector("list", length(timepoints))
names(results.lmm) <- timepoints

for (tp in timepoints) {
  results.lm[[tp]] <- read.table(file=paste0("QTL/LM_", tp, ".txt"), sep = '\t')
  results.lmm[[tp]] <- read.table(file=paste0("QTL/LMM_", tp, ".txt"), sep = '\t')
}

plot(c(0, 150), c(0, 30), t = 'n')
for (tp in timepoints) {
 points(-log10(results.lm[[tp]][1:150,"MF"]), t ='l')
}
points(-log10(LMMglobal["F_M_C",]), t = 'l', col='blue', lwd=2)
points(-log10(LMMglobal["F_M",]), t = 'l', col='red', lwd=2)

#write.table(t(LMMglobal), file=paste0("QTL/LMM_Global.txt"), sep='\t', quote=FALSE)

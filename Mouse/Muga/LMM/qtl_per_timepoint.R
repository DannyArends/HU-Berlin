# QTL mapping per timepoint of the MegaMuga data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
library(parallel)

setwd("D:/Edrive/Mouse/DNA/MegaMuga/inputF2")
genotypes <- read.table(file="cleaned_recoded_genotypes_F2.txt", sep = '\t', check.names = FALSE)
map <- read.table(file="cleaned_map_25MbGap.txt", sep = '\t')
phenotypes <- read.table(file="cleaned_phenotypes_F2.txt", sep = '\t')

# Extract covariates
timepoints <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")
subfamily <- as.factor(phenotypes[, "Vater"]) # Fixed effect: Subfamily structure (factor)
littersize <- as.numeric(phenotypes[, "WG2"]) # Fixed effect: Size of the litter (linear effect)
litternumber <- as.factor(as.numeric(phenotypes[, "W.Label"] == "A") + 1)
ltype <- factor(paste0(litternumber, "_", littersize), levels = c("1_8", "1_9", "1_12", "2_8", "2_9", "2_10", "2_11", "2_12"))

season <- as.factor(phenotypes[, "Season"]) # Fixed effect: Season when born (factor)
topmarker <- as.factor(unlist(genotypes["UNC5048297", ])) # topmarker 

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

#x <- unlist(genotypes["UNC31553261", ])
#MF <- factor(as.character(x), levels = c("B", "H", "N"))
#phenotype = as.numeric(phenotypes[, "d21"])

### Map everything using a per timepoint linear model
results.lm <- vector("list", length(timepoints))
names(results.lm) <- timepoints
for(tp in timepoints){
  cl <- makeCluster(getOption("cl.cores", 6))
  clusterExport(cl, "subfamily") # export subfamily to the nodes
  clusterExport(cl, "ltype") # export littersize to the nodes
  clusterExport(cl, "topmarker") # export season to the nodes
  res <- parApply(cl, genotypes, 1, function(x, phenotype){
    MF <- factor(as.character(x), levels = c("N", "H", "B"))
    model.lm <- lm(phenotype ~ subfamily + ltype + MF)
    model.lmc <- lm(phenotype ~ subfamily + ltype + topmarker + MF)
    pvals <- c(rep(NA, 4), rep(0, 6))
    names(pvals) <- c("subfamily", "ltype", "MF", "MF_C", "MFB", "MFH", "MFN", "MFCB", "MFCH", "MFCN")
    pvals[1:3] <- anova(model.lm)[1:3, "Pr(>F)"]
    pvals[4] <- anova(model.lmc)[4, "Pr(>F)"]
    pvals[5:7] <- model.lm$coefficients[c("MFB", "MFH", "MFN")]
    pvals[8:10] <- model.lmc$coefficients[c("MFB", "MFH", "MFN")]
    pvals[names(which(is.na(pvals[5:10])))] <- 0
    mG <- names(which(table(MF) == 0))
    if(length(mG) > 0){
      pvals[paste0("MF",mG)] <- NA
      pvals[paste0("MFC",mG)] <- NA
    }
    return(pvals)
  }, phenotype = as.numeric(phenotypes[, tp]))
  stopCluster(cl)
  rownames(res) <- c("subfamily", "ltype", "MF", "MF_C", "MFB", "MFH", "MFN", "MFCB", "MFCH", "MFCN")
  write.table(t(res), file=paste0("QTL/LM_", tp, ".txt"), sep='\t', quote=FALSE)
  results.lm[[tp]] <- res
  cat("Done ", tp, "\n")
}

#plot(-log10(results.lm[[1]]["MF",]), col=c("orange", "blue", "green")[apply(results.lm[[1]][c("MFB", "MFH", "MFN"),],2,which.max)], pch=19)

plot(c(0, max(map[, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n',las=2)
for (tp in timepoints) {
  for (chr in names(chrs)) {
    onChr <- which(map[,"Chr"] == chr)
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]]["MF", onChr]), t ='l', col='gray',lwd=3)
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]]["MF_C", onChr]), t ='l', col='blue',lwd=2)
  }
}
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")

### Map everything using a per timepoint linear mixed model, using subfamily as a random effect
results.lmm <- vector("list", length(timepoints))
names(results.lmm) <- timepoints
for(tp in timepoints){
  cl <- makeCluster(getOption("cl.cores", 6))
  clusterExport(cl, "subfamily") # export subfamily to the nodes
  clusterExport(cl, "ltype") # export littersize to the nodes
  clusterExport(cl, "topmarker") # export season to the nodes
  clusterEvalQ(cl, library(lme4)) # load lme4 per node
  res <- parApply(cl, genotypes, 1, function(x, phenotype){
    ctrl = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun=20000))
    noG <- which(is.na(x))
    MF <- factor(as.character(x), levels = c("N", "H", "B"))
    if(length(noG) > 0){
      phenotype <- phenotype[-noG]
      MF <- MF[-noG]
      subfamily <- subfamily[-noG]
      ltype <- ltype[-noG]
      topmarker <- topmarker[-noG]
    }
    model.full <- lmer(phenotype ~ (1|subfamily) + ltype + MF, control=ctrl, REML = FALSE)
    model.null <- lmer(phenotype ~ (1|subfamily) + ltype, control=ctrl, REML = FALSE)
    model.fullC <- lmer(phenotype ~ (1|subfamily) + ltype + topmarker + MF, control=ctrl, REML = FALSE)
    model.nullC <- lmer(phenotype ~ (1|subfamily) + ltype + topmarker, control=ctrl, REML = FALSE)
    pvals <- c(rep(NA, 2), rep(0, 6))
    names(pvals) <- c("Full", "Full_C", "MFB", "MFH", "MFN", "MFCB", "MFCH", "MFCN")
    pvals[1] <- anova(model.full, model.null)["model.full", "Pr(>Chisq)"]
    pvals[2] <- anova(model.fullC, model.nullC)["model.fullC", "Pr(>Chisq)"]
    pvals[3:5] <- summary(model.full)$coefficients[, "Estimate"][c("MFB", "MFH", "MFN")]
    pvals[6:8] <- summary(model.fullC)$coefficients[, "Estimate"][c("MFB", "MFH", "MFN")]
    pvals[names(which(is.na(pvals[3:8])))] <- 0
    mG <- names(which(table(MF) == 0))
    if(length(mG) > 0){
      pvals[paste0("MF",mG)] <- NA
      pvals[paste0("MFC",mG)] <- NA
    }
    return(pvals)
  }, phenotype = as.numeric(phenotypes[, tp]))
  stopCluster(cl)
  rownames(res) <- c("Full", "Full_C", "MFB", "MFH", "MFN", "MFCB", "MFCH", "MFCN")
  write.table(t(res), file=paste0("QTL/LMM_", tp, ".txt"), sep='\t', quote=FALSE)
  results.lmm[[tp]] <- res
  cat("Done ", tp, "\n")
}

tpcolors <- colorRampPalette(c("lightblue", "darkblue"))(length(timepoints))
tpcolorsC <- colorRampPalette(c("firebrick1", "firebrick4"))(length(timepoints))
names(tpcolors) <- timepoints
names(tpcolorsC) <- timepoints

plot(c(0, max(map[, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n',las=2)
for (tp in timepoints) {
  for (chr in names(chrs)) {
    onChr <- which(map[,"Chr"] == chr)
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][onChr, "MF"]), t ='l', col = tpcolors[tp])
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][onChr, "MF_C"]), t ='l', col = tpcolorsC[tp])
  }
}
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")


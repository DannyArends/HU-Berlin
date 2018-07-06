setwd("D:/Edrive/Horse/DNA/Analysis2017")

### Write out the filtered genotype data
map <- read.table("genotypes_map_filtered.txt", sep="\t",header=TRUE)
genotypes <- read.table("genotypes_snp_filtered.txt", sep="\t",na.strings=c("", "NA"))
phenotypes <- read.table("phenotypes_corrected_filtered.txt", sep="\t")
pheALL <- read.table("phenotypes.txt", sep="\t")

## find outliers and put them to NA
for(trait in rownames(phenotypes)) {
  outlierLOW <- which(as.numeric(phenotypes[trait,]) < mean(as.numeric(phenotypes[trait,])) - 2*sd(as.numeric(phenotypes[trait,])))
  outlierHIGH <- which(as.numeric(phenotypes[trait,]) > mean(as.numeric(phenotypes[trait,])) + 2*sd(as.numeric(phenotypes[trait,])))
  outliers <- c(outlierLOW, outlierHIGH)
  if(length(outliers) > 0) phenotypes[trait, outliers] <- NA
}

strains <- as.factor(unlist(pheALL["Strain",colnames(phenotypes)]))

phenotypes.corrected <- phenotypes

mtable <- cbind(apply(phenotypes[,names(strains[which(strains == "S")])],1,median,na.rm=TRUE),apply(phenotypes[,names(strains[which(strains == "S")])],1,sd,na.rm=TRUE),
  apply(phenotypes[,names(strains[which(strains == "K")])],1,median,na.rm=TRUE), apply(phenotypes[,names(strains[which(strains == "K")])],1,sd,na.rm=TRUE),
  apply(phenotypes[,names(strains[which(strains == "H")])],1,median,na.rm=TRUE), apply(phenotypes[,names(strains[which(strains == "H")])],1,sd,na.rm=TRUE))

colnames(mtable) <- c("S(mean)", "S(SD)", "K(mean)", "K(SD)", "H(mean)", "H(SD)")
mtable <- cbind(mtable, "Pval" =  NA)

for(trait in rownames(phenotypes)) {
  model <- lm(as.numeric(phenotypes[trait, ]) ~ strains)
  mtable[trait, "Pval"] <- round(anova(model)[[5]][1],3)
}

write.table(mtable, "Table_withSD.txt", sep="\t", quote=FALSE)

for(trait in rownames(phenotypes)) {
  pval <- anova(lm(as.numeric(phenotypes[trait, ]) ~ strains))[[5]][1]
  cat(trait, " ", pval, "\n")
  if(pval < 0.1){
    cat("Correcting", trait, "\n")
    model <- lm(as.numeric(phenotypes[trait, ]) ~ strains)
    newpheno <- model$residuals + model$coefficients["(Intercept)"]
    phenotypes.corrected[trait, as.numeric(names(newpheno))] <- round(newpheno, d = 2)
  }
}

phenotypes <- phenotypes.corrected 

write.table(phenotypes.corrected, "phenotypes_filtered_NoSex_Outliers_Strain.txt", sep="\t", quote=FALSE)

# Do the GWAS using a single QTL model *Do not correct for race of horse
pvalues <- matrix(NA, nrow(genotypes), nrow(phenotypes), dimnames = list(rownames(genotypes), rownames(phenotypes)))
for(trait in rownames(phenotypes)) {
  phe <- as.numeric(phenotypes[trait,])
  for(marker in rownames(genotypes)){
    pvalues[marker, trait] <- anova(lm(phe ~ as.factor(as.character(genotypes[marker, ]))))[[5]][1]
  }
  cat("Done", trait, "\n")
}

allelefreq <- apply(genotypes, 1 , function(x){
  tbl <- table(unlist(lapply(x, strsplit, "")))
  min(tbl / sum(tbl))
})

write.table(pvalues, "pvalues.txt", sep="\t", quote=FALSE)
write.table(round(-log10(pvalues),d=2), "lodscores.txt", sep="\t", quote=FALSE)

pvalues <- read.table("pvalues.txt", sep = "\t")

chr <- as.numeric(as.character(map$Chromosome))
chr[is.na(chr)] <- 32

padjusted <- apply(pvalues,2,function(x){
  p.adjust(x, "fdr")
})

padjustedBF <- apply(pvalues,2,function(x){
  p.adjust(x, "bonferroni")
})

lodscores <- -log10(pvalues)
lodadjusted <- -log10(padjusted)

for(trait in rownames(phenotypes)) {
  significant <- which(padjusted[,trait] < 0.1)
  if(length(significant) > 0){
    cat(trait, " ", length(significant), "\n")
  }
  cat(trait, " ", max(lodscores[,trait]), " ", min(pvalues[,trait]), " ", min(padjusted[,trait]), "\n")
}
#op <- par(mfrow=c(3,3))

for(trait in rownames(phenotypes)) {
  # NG, ChW
  trait <- "ChW"
  significant <- which(padjusted[,trait] < 0.15)

  cbind(map[significant,], MAF = round(map[significant,"MAF"], d = 2), P = pvalues[significant, trait], "P(FDR)" = padjusted[significant, trait] , "P(B)" = padjustedBF[significant, trait])

  observed.p <- pvalues[, trait]
  expected.p <- (rank(observed.p, ties.method="first")+0.5) / (length(observed.p) + 1)

  expected.l <- -log10(sort(expected.p))
  observed.l <- -log10(sort(observed.p))

  plot(x = expected.l, y = observed.l, xlim=c(0,max(c(expected.l, observed.l))), ylim=c(0,max(c(expected.l, observed.l))))
  abline(0, 1)

  # Estimate the lambda value
  lambda <- median(qchisq(1.0 - observed.p, 1)) /  qchisq(0.5, 1)
  cat(trait, " ", lambda, "\n")
}
trait <- "NG"
layout(matrix(c(1,1,2,3),2,2,byrow=T))
chrl <- c()
for(x in 1:32){
  chrl <- c(chrl, max(map[which(chr == x),"Position"]) + 15000000)
}
op <- par(mar=c(0, 4, 4, 2) + 0.1)
chrss <- c(0, chrl)
#for(trait in rownames(phenotypes)) {
  plot(y=c(0,6), x=c(-10000000, sum(as.numeric(chrl))-5000000), t = 'n', xlab="", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main="Manhattan plot - Neck Girth")
  mm <- c()
  for(x in 1:32){
    idx <- which(chr == x)
    points(x = sum(chrss[1:x]) + as.numeric(map[idx, "Position"]), y = lodscores[idx,trait], t = 'h', col=c("black", "cornflowerblue")[(x %% 2 == 1) + 1])
    m <- sum(chrss[1:x]) + ((chrl[x]-15000000) / 2)
    if(x < 32)text(x= m, y = 0.20, x, col=c("white", "black")[(x %% 2 == 1) + 1],cex=0.8)
    if(x == 32) text(x= m, y = 0.20, "X", col=c("white", "black")[(x %% 2 == 1) + 1],cex=0.8)
    cat(chrss[x], lc, m, "\n")
    l <- l + lc
    mm <- c(mm,m)
  }
  abline(h = 5, col = "orange", lty=2)
  abline(h = 5.5, col = "green", lty=2)
  #axis(1, at=mm, c(seq(1,31,1),"X"))
#}
  observed.p <- pvalues[, trait]
  expected.p <- (rank(observed.p, ties.method="first")+0.5) / (length(observed.p) + 1)

  expected.l <- -log10(sort(expected.p))
  observed.l <- -log10(sort(observed.p))

  # Estimate the lambda value
  lambda <- median(qchisq(1.0 - observed.p, 1)) /  qchisq(0.5, 1)
  cat(trait, " ", lambda, "\n")
  op <- par(mar=c(5, 4, 4, 2) + 0.1)
  plot(x = c(0,5), y = c(0,6), xlab="Expected -log10(p)", ylab="Observed -log10(p)", main="Neck Girth", xaxs="i", yaxs="i", yaxt='n')
  points(x = expected.l, y = observed.l, pch=19, col=c("orange", rep("cornflowerblue",(length(expected.l) - 1))))
  mtext(paste0("Genomic inflation factor λ=", round(lambda, 3)),cex=0.8)
  axis(2, at=c(0:6), c(0:6), las=2)
  abline(0, 1)

  boxplot(as.numeric(phenotypes[trait,]) ~ as.character(unlist(genotypes["BIEC2_772752",])), main="Neck Girth @ BIEC2_772752", yaxt='n',xlab="Genotype", ylab="Neck Girth (cm)")
  axis(2, at=c(105,110,115,120), c(105,110,115,120), las=2)
  mtext(paste0("Minor allele T MAF=", round(map["BIEC2_772752", "MAF"], 2)), cex=0.8)
  
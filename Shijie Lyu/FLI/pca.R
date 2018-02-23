library(RColorBrewer)
library(biomaRt)
library(parallel)

setwd("D:/Edrive/Chicken/ClassicalPhenotypes/FLIData")

mdata <- readLines("Asia_Game_vs_Bantam_chicken.txt")

alldata <- strsplit(mdata, "\t")
dmatrix <- matrix(unlist(alldata), nrow=length(alldata), ncol=length(alldata[[1]]), byrow=TRUE)
dmatrix <- t(dmatrix)

colheader <- dmatrix[1,][-1]
rowheader <- dmatrix[,1][-1]
dmatrix <- dmatrix[-1,-1]

colnames(dmatrix) <- colheader
rownames(dmatrix) <- rowheader

dmatrix[dmatrix == "?_?"] <- NA
dmatrix[dmatrix == "?"] <- NA

game <- c("ASrb", "MAxx", "OFrbx", "IKxx")#, "SHsch")                            # more breeds: SHsch
bantam <- c("CHgesch", "CHschw", "KSgw", "OHgh", "OHsh")#, "ZCsch", "ZCw")       # more breeds: ZCsch, ZCw

# Colors we'll use throughout all the plots
mycolors <- brewer.pal(length(c(game,bantam)),"Paired")
names(mycolors) <- c(game,bantam)

markers <- rownames(dmatrix)[which(grepl("AX-", rownames(dmatrix)))]
individuals <- colnames(dmatrix)[which(dmatrix["popPk", ] %in% c(game,bantam))]

# Phenotypes we are allowed to use
phenonames <- c("popPk", "sex10", "AlterMo", "gew", "armL", "armR", "laufL", "laufR", "laufdiL", "laufdiR", "kiel")

# Get the genotypes and phenotypes in their own matrices
genotypes <- dmatrix[markers, individuals]
phenotypes <- dmatrix[phenonames, individuals]

chromosomes <- 1:29
names(chromosomes) <- c(1:28, "Z")

# Load the marker annotation
amatrix <- read.csv("GalGal5_Claudia.txt", row.names=1, colClasses="character")
annot <- amatrix[which(amatrix[,"Chromosome"] %in% names(chromosomes)),c("Chromosome", "Position")]

# Call rate: Keep SNPs with a callrate lower or equal to 5%
callrate <- apply(genotypes, 1, function(x){ return(sum(is.na(x)) / length(x)) })
genotypes <- genotypes[which(callrate <= 0.05), ]

# SNP minor allele frequencies
mafs <- apply(genotypes,1,function(x){
  mt <- table(unlist(strsplit(x, "_")))
  min(mt / sum(mt))
})

# MAFs above (or equal) to 5% MAF & not monomorphic (MAF = 1)
genotypes <- genotypes[which(mafs >= 0.05 & mafs <= 0.5), ]

# Make sure the markers have annotation
genotypes <- genotypes[which(rownames(genotypes) %in% rownames(annot)),]
annot <- annot[which(rownames(annot) %in% rownames(genotypes)),]

# Order the markers per chromosome
genotypes <- genotypes[rownames(annot), ]

numGT <- genotypes

# Recode to numeric, so we can do distances
numGT[numGT == "A_A"] <- 1
numGT[numGT == "A_B"] <- 2
numGT[numGT == "B_B"] <- 3

numGT <- apply(numGT, 2, function(x){as.numeric(x)} )
rownames(numGT) <- rownames(genotypes)

ind.distances <- dist(t(numGT[sample(rownames(numGT), 100000), ]))
dendro <- as.dendrogram(hclust(ind.distances, method="complete"))

# function to color and rename labels in the dendrogram
colLabels <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    breed <- dmatrix["popPk", a$label]
    mcol = mycolors[breed]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = as.character(mcol))
    attr(n, "label") <- paste0(breed) #,"_", a$label)
  }
  n
}

# Color the dendrogram
colDendro = dendrapply(dendro, colLabels)

# Plot the dendrogram
png(paste0("Dendrogram.png"), width=2048, height=1024, res = 300, pointsize = 6)
  par(cex=0.5)
  plot(colDendro)
dev.off()

# Found outliers, remove them from genotypes, numGT and phenotypes
outliers <- c(1252, 1331, 1347, 2136, 2214)

markersOnZ <- rownames(annot)[which(annot[,"Chromosome"] == "Z")]
gtTables <- apply(numGT[markersOnZ,],2, function(x){
  table(x)
})

shouldBeFemale <- names(gtTables)[which(phenotypes["sex10", names(gtTables)] == "0")] # 1 = Male, 0 = Female
shouldBeMale <- names(gtTables)[which(phenotypes["sex10", names(gtTables)] == "1")] # 1 = Male, 0 = Female

names(which(is.na(lapply(gtTables[shouldBeMale], "[", "2"))))
names(which(!is.na(lapply(gtTables[shouldBeFemale], "[", "2"))))

# Female which is actually a male: 1141, 1252 (outlier), 2152
# Male which is actually a Female: 1331 (outlier)
phenotypes["sex10", "1141"] <- 1
phenotypes["sex10", "1252"] <- 1
phenotypes["sex10", "2152"] <- 1
phenotypes["sex10", "1331"] <- 0

# Remove the outlier as seen in the dendrogram
numGT <- numGT[, which(!colnames(numGT) %in% outliers)]
genotypes <- genotypes[, which(!colnames(genotypes) %in% outliers)]
phenotypes <- phenotypes[, which(!colnames(phenotypes) %in% outliers)]

numGTblank <- apply(numGT, 1, function(x){
  tbl <- table(x)
  tooSmall <- names(which(tbl <= 10))
  if(length(tooSmall) > 0){
    x[x %in% tooSmall] <- NA
  }
  return(x)
})
numGT <- t(numGTblank)
genotypes[is.na(numGT)] <- NA  # Graft the blanks into the character genotypes

# Recalculate SNP minor allele frequencies
mafs <- apply(genotypes, 1, function(x){
  mt <- table(unlist(strsplit(x, "_")))
  min(mt / sum(mt))
})

# Reapply our MAFs filter
genotypes <- genotypes[names(which(mafs >= 0.05 & mafs <= 0.5)), ]
numGT <- numGT[names(which(mafs >= 0.05 & mafs <= 0.5)), ]

# Recalculate the callrate
callrate <- apply(numGT, 1, function(x){ sum(is.na(x)) / length(x) })

# Transpose the numGT matrix
numGT <- t(numGT)
# genotypes for PCA (no missing data)
numGTforPCA <- numGT[,which(callrate==0)]

# Which individuals are which type
ind.bantam <- colnames(dmatrix)[which(dmatrix["popPk", ] %in% bantam)]
ind.bantam <- ind.bantam[which(ind.bantam %in% rownames(numGT))]

ind.game <- colnames(dmatrix)[which(dmatrix["popPk", ] %in% game)]
ind.game <- ind.game[which(ind.game %in% rownames(numGT))]

### Scale down the phenotypes and QTL mapping
phenotypes <- phenotypes[, rownames(numGT)]

pheno <- phenotypes[c("gew", "kiel"), ]

pheno <- rbind(pheno, arm = (as.numeric(phenotypes["armL", ]) + as.numeric(phenotypes["armR", ])) / 2)
pheno <- rbind(pheno, lauf = (as.numeric(phenotypes["laufL", ]) + as.numeric(phenotypes["laufR", ])) / 2)
pheno <- rbind(pheno, laufdi = (as.numeric(phenotypes["laufdiL", ]) + as.numeric(phenotypes["laufdiR", ])) / 2)

population <- as.factor(phenotypes["popPk",])
sex <- as.factor(phenotypes["sex10",])
age <- as.numeric(phenotypes["AlterMo",])

phenoNames <- rownames(pheno)

# Adjust phenotypes for the (corrected) sex of the animal, within each population
pheno.adj <- matrix(NA, nrow(pheno), ncol(pheno), dimnames = list(rownames(pheno), colnames(pheno)))
adj.LOD <- matrix(NA, length(unique(population)), length(phenoNames), dimnames = list(unique(population), phenoNames))
for(pop in unique(population)){
  for(x in phenoNames) {
    inPop <- names(population)[which(population == pop)]
    model <- lm(as.numeric(pheno[x,inPop]) ~ sex[inPop])
    pheno.adj[x, inPop] <- round(model$coefficients["(Intercept)"] + model$residuals,2)
    cat(pop, " ", x, " ", -log10(anova(model)[[5]][1]), "\n")
    adj.LOD[pop, x] <- -log10(anova(model)[[5]][1])
  }
  cat("DONE", pop, "\n")
}

# Calculate the cummultative position of markers for the plot
chr.lengths <- c()
chr.starts <- c(0)
for(x in 1:length(chromosomes)){
  chrlength <- max(as.numeric(annot[annot[, "Chromosome"] == names(chromosomes)[x],"Position"]))
  chr.lengths <- c(chr.lengths, chrlength)
  chr.starts <- c(chr.starts, chr.starts[x] + chrlength + 5000000)
  cat(x," ", names(chromosomes)[x], " ", chrlength, "\n")
}

marker.chr <- as.numeric(chromosomes[annot[colnames(numGT), "Chromosome"]])
marker.pos <- as.numeric(annot[colnames(numGT), "Position"]) + chr.starts[marker.chr]
chr.col = c("orange", "gray")[(marker.chr %% 2) + 1]

# QTL mapping, with population structure correction
for(phe in rownames(pheno.adj)){
  # Map the QTLs for this phenotype
  cl <- makeCluster(4)
  clusterExport(cl, "pheno.adj")
  clusterExport(cl, "population")
  qtl <- parApply(cl, numGT, 2, function(geno, phe = "gew"){
    return(na.omit(anova(lm(pheno.adj[phe,] ~ as.factor(population) + as.factor(geno)))[[5]]))
  }, phe = phe)
  stopCluster(cl)

  qtlTable <- cbind(population = unlist(lapply(qtl, "[", 1)), additive = unlist(lapply(qtl, "[", 2)))

  write.table(qtlTable, file=paste0("all_qtl_", phe, ".txt"), sep = "\t", quote = FALSE)

  # Adjust LOD scores using BH correction
  lodscores.adj <- -log10(p.adjust(qtlTable[,"additive"], "BH"))

  png(paste0("plot_qtl_", phe, ".png"), width=1024, height=600)
    plot(x = marker.pos, y = lodscores.adj, col=chr.col, pch=19, cex=0.7, xaxt='n', las=2, xlab="Chromosome", ylab="LOD adjusted(BH)", main=paste0("Phenotype: ",phe))
    axis(1, at = chr.starts[chromosomes] + (chr.lengths /2), chromosomes)
    abline(h=-log10(0.10), col="orange")
    abline(h=-log10(0.05), col="green", lty=2)
    abline(h=-log10(0.01), col="green")
  dev.off()

  # Significant SNPs after BH correction
  significant.snps <- which(lodscores.adj > -log10(0.05))
  cat(phe, " found ", length(significant.snps), " significant SNPs\n")
  
  png(paste0("plot_QQ_", phe, ".png"), width=800, height=800)
  # Inflation factor
  rankedData <- rank(unlist(qtlTable[,"additive"]), ties.method="first") + 0.5
  plot(-log10(rankedData / (max(rankedData) + 1)), -log10(unlist(qtlTable[,"additive"])))
  abline(a=0, b=1)
  lambda <- round(median(qchisq(1.0 - unlist(qtlTable[,"additive"]), 1),na.rm=TRUE) /  qchisq(0.5, 1),3)
  cat(phe, " Lambda:", lambda, "\n")
  # Pvalues to Zscores
  zscores <- qnorm(1 - (qtlTable[,"additive"] / 2)) 
  
  # Zscored deflated for any inflation in lamdba
  p.deflated = pchisq((zscores^2)/lambda, df = 1, lower = FALSE)
  
  lambda <- round(median(qchisq(1.0 - p.deflated, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)
  # Now lambda should be exactly 1 !
  cat(phe, " Lambda deflated:", lambda, "\n")
  rankedData <- rank(p.deflated, ties.method="first") + 0.5
  points(-log10(rankedData / (max(rankedData) + 1)), -log10(p.deflated), col="blue")
  dev.off()
  p.deflated.adj = p.adjust(p.deflated, "BH")
  
  significant.snps.deflated <- which(-log10(p.deflated.adj) > -log10(0.05))
  cat(phe, " found ", length(significant.snps.deflated), " significant SNPs after deflation\n")
  sign.table <- cbind(annot[names(significant.snps),],
                      maf = round(100 * mafs[names(significant.snps)], 1), 
                      lod = round(-log10(qtlTable[significant.snps, 2]),2), 
                      lod.BH.adj = round(lodscores.adj[significant.snps],2),
                      lod.deflated = round(-log10(p.deflated[significant.snps]),2), 
                      lod.deflated.BH.adj = round(-log10(p.deflated.adj[significant.snps]),2)                      
                     )
  write.table(sign.table, file=paste0("sign_qtl_", phe, ".txt"), sep = "\t", quote = FALSE) 
  
}

qtlALL <- NULL
for(phe in rownames(pheno.adj)){
  qtlTable <- read.table(file=paste0("all_qtl_", phe, ".txt"), sep = "\t")
  qtlALL <- cbind(qtlALL, qtlTable[,2])
}
colnames(qtlALL) <- rownames(pheno.adj)
rownames(qtlALL) <- rownames(qtlTable)
qtlALL <- -log10(qtlALL)

# PCA analysis
pca <- prcomp(numGTforPCA)

breed <- dmatrix["popPk", names(pca$x[,1])]
mcol = mycolors[breed]
type <- c(15, 18)[as.numeric(breed %in% game) + 1]
png(paste0("plot_PC1_PC2.png"), width=1024, height=800, res=450, pointsize=6)
  par(cex=0.4)
  plot(pca$x[,1], pca$x[,2], col=mcol, pch=type, xlab="Principal component 1", ylab="Principal component 2")
  legend("topright", names(mycolors), col=mycolors, pch=c(15, 18)[as.numeric(names(mycolors) %in% game) + 1])
dev.off()

# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){ return(var.loadings*comp.sdev) }

# Contribution of the variables to the principal component
contrib <- function(var.cos2, comp.cos2){ return((var.cos2 * 100)/comp.cos2) }

# Variable correlation/coordinates
var.coord <- t(apply(pca$rotation, 1, var_cor_func, pca$sdev))
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2, 1, contrib, comp.cos2))

# Pvalue and LOD scores for contribution
pcaDistObs <- (sign(pca$rotation[,1]) * var.contrib[,1])    #  mydistExp <- (sign(pca$rotation[,1]) * sqrt(var.contrib[,1]))
pcaPvals <- 2 * pnorm(abs(pcaDistObs), 0, sd=sd(pcaDistObs), lower.tail=FALSE)
pcaPvals.adj <- p.adjust(pcaPvals, "holm")

inPCA <- rownames(annot[rownames(var.contrib),])
#inPCAgood <- names(which(mafs[inPCA] > 0.3))
annotPCA <- annot[inPCA,]


plot(c(0, 180000), c(-10,10), t = 'n')
points(-log10(pcaPvals))
points(-runmed(-log10(pcaPvals), 15))


# Significance
significantPCA <- names(which(pcaPvals.adj <= 0.1))
length(significantPCA)


plot(c(0,30), c(0, max(chr.lengths)), t ='n',yaxt='n', xaxt='n', ylab="Position (Mbp)", xlab="Chromosome")
axis(1, at=1:length(chromosomes), names(chromosomes))
axis(2, at=seq(0,max(chr.lengths), 25000000), seq(0,max(chr.lengths), 25000000)/1000000, las=2)
cnt <- 1
aa <- lapply(chr.lengths, function(x){ lines(c(cnt,cnt),c(0, chr.lengths[cnt])); cnt <<- cnt+1 })
aa <- apply(annot[significantPCA,],1,function(x){ points(chromosomes[x[1]],x[2], pch="-", col='orange', cex=2)})

pca.table <- cbind(annot[significantPCA,], 
                   contribution = round(var.contrib[significantPCA,1],4), 
                   MAF = round(mafs[significantPCA],2), 
                   LOD = round(-log10(pcaPvals)[significantPCA], 1))

write.table(pca.table, file=paste0("pca_sign.txt"), sep = "\t", quote = FALSE) 


plot(mafs[rownames(var.contrib)], round(pcaLOD[rownames(var.contrib)], 1))
regions.size <- 1000000
results <- NULL
for(chr in names(chromosomes)){
  chr.length <- chr.lengths[chromosomes[chr]]
  region.start <- 0
  onChr <- rownames(annotPCA)[which(annotPCA[,"Chromosome"] == chr)]
  while(region.start < (chr.length - regions.size)){
    inRegion <- rownames(annotPCA[onChr,])[which(as.numeric(annotPCA[onChr, "Position"]) >= region.start & as.numeric(annotPCA[onChr, "Position"]) <= (region.start + regions.size))]
    results <- rbind(results, c(chr, region.start+(regions.size/2),length(inRegion), mean(pcaPvals[inRegion])))
    region.start <- region.start + regions.size
  }
  cat("Done chromosome ", chr, "\n")
}


GHR <- rownames(annotPCA)[which(annotPCA[,"Chromosome"] == "Z" & as.numeric(annotPCA[,"Position"]) >= 13000000 & as.numeric(annotPCA[,"Position"]) <= 14000000)]


# QTL mapping body weight, without population structure correction
cl <- makeCluster(4)
clusterExport(cl, "pheno.adj")

qtl <- parApply(cl, numGTforPCA, 2, function(geno, phe = "gew"){
  model <- lm(pheno.adj[phe,] ~ geno)
  return(sign(model$coefficients["geno"]) * na.omit(anova(model)[[5]]))
}, phe = "gew")
stopCluster(cl)

# Convert to LOD scores
lodsQTL <- -log10(abs(qtl))

mydistObs <- (sign(pca$rotation[,1]) * var.contrib[,1])    #  mydistExp <- (sign(pca$rotation[,1]) * sqrt(var.contrib[,1]))
myps <- 2 * pnorm(abs(mydistObs), 0, sd=sd(mydistObs), lower.tail=FALSE)
plot(-log10(myps), col=as.numeric(rownames(var.contrib) %in% snpsHigh)+1 , pch=19, cex=0.5)
abline(h = -log10(0.05 / nrow(var.contrib)))

# Convert to LOD scores
lodsPCA <- (-log10(myps))
model <- lm(lodsQTL ~ lodsPCA)
QTLcutoff <- predict(model, data.frame(lodsPCA=-log10(0.05 / nrow(var.contrib))))
plot(x = lodsPCA, y = lodsQTL, cex=0.5, pch=19)
abline(a=model$coefficients["(Intercept)"], b = model$coefficients["lodsPCA"], col='green')
abline(h = QTLcutoff, col='red')
abline(v = -log10(0.05 / nrow(var.contrib)), col='blue')

regionsQTL <- read.table("regionsQTL.txt", sep="\t", header=TRUE)
regions <- paste0(regionsQTL[,1],":",regionsQTL[,2],":",regionsQTL[,4])
biomart <- useMart("ENSEMBL_MART_ENSEMBL", "ggallus_gene_ensembl")

for(r in regions){
  res.biomart <- getBM(c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position", "strand", "description"), 
                      filters="chromosomal_region", values=r, mart=biomart)
  res.biomart <- res.biomart[which( res.biomart[,2] != ""),]
  res.biomart[,7] <- unlist(lapply(strsplit(res.biomart[,7], " [",fixed=TRUE),"[",1))
  write.table(res.biomart, file=paste0(gsub(":", "_",r),".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}

regionsPCA <- read.table("regionsPCA.txt", sep="\t", header=TRUE)
regions <- paste0(regionsPCA[,1],":",regionsPCA[,2],":",regionsPCA[,4])

for(r in regions){
  res.biomart <- getBM(c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position", "strand", "description"), 
                      filters="chromosomal_region", values=r, mart=biomart)
  res.biomart <- res.biomart[which( res.biomart[,2] != ""),]
  res.biomart[,7] <- unlist(lapply(strsplit(res.biomart[,7], " [",fixed=TRUE),"[",1))
  write.table(res.biomart, file=paste0("PCA_", gsub(":", "_",r),".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}

ghrLoc <- c(13478526, 13558649)
ghrSNPs <- rownames(annotPCA[annotPCA$Chromosome == "Z" & as.numeric(annotPCA$Position) > ghrLoc[1] & as.numeric(annotPCA$Position) < ghrLoc[2],])

apply(numGTforPCA[bantamIND,ghrSNPs],2,table)
apply(numGTforPCA[gameIND,ghrSNPs],2,table)

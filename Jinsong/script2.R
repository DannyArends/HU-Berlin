# Script by Danny Arends, for analysis of BXD bone density data by glm

library("lme4")

setwd("D:/Ddrive/Collegues/Jinsong")

# Genotypes
genotypes <- read.table("BXD.geno", sep="\t", header=TRUE, row.names=2, colClasses="character", na.string = c("U", "H", "NA"))
genotypes <- cbind(genotypes, "DBA/2J" = rep("D", nrow(genotypes)))
genotypes <- cbind(genotypes, "C57BL/6J" = rep("B", nrow(genotypes)))

# Genetic map
map <- genotypes[,1:3]

# Length of chromosomes, and cumulative positions
chrs <- as.character(c(seq(1, 19, 1), "X"))
chrs.lengths <- NULL
chrs.starts <- 0
for(i in 1:length(chrs)) {
  chrs.lengths <- c(chrs.lengths, max(as.numeric(map[which(map[,"Chr"] == chrs[i]), "cM"]), na.rm = TRUE))
  chrs.starts <- c(chrs.starts, 100 + chrs.starts[i] + chrs.lengths[i])
}
names(chrs.starts) <- chrs
map <- cbind(map, cMplot = (as.numeric(map[,"cM"]) + chrs.starts[map[,"Chr"]]))

# Phenotypes
phenotypes <- read.table("BXD.bone", sep="\t", header=TRUE, colClasses="character", check.names=FALSE)

genotypes[1:5,1:10]; dim(genotypes)
phenotypes[1:5,1:10]; dim(phenotypes)

# Phenotype names of interest, the rest is covar
phe.names <- colnames(phenotypes[-c(1:5)])

# Helper functions, simple linear models to get a 'handle; on the data
get.eff <- function(phenotypes, phe.name = "Femur Length (mm) F", effect = "Sex") {
  phenotype.data <- as.numeric(phenotypes[, phe.name])
  group.factor <- as.factor(phenotypes[, effect])
  boxplot(phenotype.data ~ group.factor, main = phe.name, sort = TRUE)
  return(lm(phenotype.data ~ group.factor))
}

get.eff.lme <- function(phenotypes, phe.name = "Femur Length (mm) F", effect = "Strain") {
  phenotype.data <- as.numeric(phenotypes[, phe.name])
  group.factor <- as.factor(phenotypes[, effect])
  boxplot(phenotype.data ~ group.factor, main = phe.name, sort = TRUE)
  return(lmer(phenotype.data ~ (1 | group.factor)))
}

get.eff.num <- function(phenotypes, phe.name = "Femur Length (mm) F", effect = "LogAGE") {
  phenotype.data <- as.numeric(phenotypes[, phe.name])
  return(lm(phenotype.data ~ as.numeric(phenotypes[, effect])))
}

get.residuals <- function(model) { 
  idx <- as.numeric(names(model$residuals)); v <- rep(NA, max(idx)); v[idx] <- as.numeric(model$residuals); return(v)
}
get.intercept <- function(model){ return(as.numeric(model$coefficients["(Intercept)"])) }

# Only take the males
males <- which(phenotypes[,"Sex"] == "F")
phe.M <- phenotypes[males, ]

# Then, adjust for the major logAge effects we observe by using the residuals + mean, on the sex adjusted data
phe.M.age.adj <- NULL
for(phe.name in phe.names){
  phe.adj <- get.residuals(get.eff.num(phe.M, phe.name, "LogAGE")) +  mean(as.numeric(phe.M[,phe.name]), na.rm = TRUE)
  phe.M.age.adj <- cbind(phe.M.age.adj, phe.adj)
}
colnames(phe.M.age.adj) <- phe.names
phe.M.age.adj <- cbind(phe.M[,1:5], phe.M.age.adj)
phe.M.age.adj[1:5, 1:8]

strainsM <- unique(phe.M[,"Strain"])
phe.M.age.adj.SA <- matrix(NA, length(strainsM), length(phe.names), dimnames=list(strainsM, phe.names))
for(strain in strainsM) {
  phe.M.age.adj.SA[strain,] <- apply(phe.M.age.adj[phe.M.age.adj[,"Strain"] == strain, phe.names], 2, function(x){
    return(mean(as.numeric(x), na.rm=TRUE))
  })
}

pheNoG <- which(!(rownames(phe.M.age.adj.SA) %in% colnames(genotypes)))
genNoP <- which(!(colnames(genotypes) %in% rownames(phe.M.age.adj.SA)))

# Keep the individuals that match between genotypes and phenotypes
phe.M.age.adj.SA <- phe.M.age.adj.SA[-pheNoG, ]
phe.M.age.adj.SA[1:5, 1:5]

genotypes <- genotypes[, - genNoP]
genotypes <- t(genotypes[,rownames(phe.M.age.adj.SA)])
genotypes[1:5,1:5]

pvalues.SA <- matrix(NA, ncol(phe.M.age.adj.SA), ncol(genotypes), dimnames=list(colnames(phe.M.age.adj.SA), colnames(genotypes)))

for(phenotype in colnames(phe.M.age.adj.SA)){
  for(marker in colnames(genotypes)) {
    pvalues.SA[phenotype, marker] <- anova(lm(phe.M.age.adj.SA[,phenotype] ~ genotypes[,marker]))[[5]][1]
  }
}

toLRS <- 2*log(10)        # Conversion constant LOD -> LRS
LRS.cutoff <- 15.00
LRS.threshold <- LRS.cutoff - (1.5 * toLRS)


marker.LRS = -log10(pvalues.SA) * toLRS

# Create plots of the different phenotypes
for(phe.name in phe.names) {
  png(paste0(gsub("%", "Pct", gsub("/", ".", gsub(" F", "", phe.name))), "_SA.png"), width = 1200, height =  600)
  plot(y = c(0, 25), x = c(0, max(map[,"cMplot"])), t = 'n', ylab = "LRS", xlab = "cM", main = paste0(gsub(" F", "", phe.name), " (females only, SA)"))
  for(i in 1:length(chrs)){
    mOnChr <- rownames(map[which(map[,"Chr"] == chrs[i]),])
    chrID <- ((i%%2) + 1)
    points(marker.LRS[phe.name, mOnChr], x=map[mOnChr,"cMplot"], t = 'h', col=c("dodgerblue", "dodgerblue4")[chrID], lwd=1)
  }
  abline(h = 15, col="green", lty = 2)
  abline(h = 12, col="gold", lty = 2)
  dev.off()
}
# Create an output table for QTL results
output.tablular <- NULL
qtl.idx <- 1
for(phe.name in phe.names) {
  for(i in 1:length(chrs)) {
    mOnChr <- rownames(map[which(map[,"Chr"] == chrs[i]),])         # Markers on chromosome
    mAboveC <- marker.LRS[phe.name, mOnChr] > LRS.cutoff            # Markers above 13.5 LRS
    if(any(mAboveC)){
      top.value <- max(marker.LRS[phe.name, mOnChr], na.rm = TRUE)     # Top LRS score on this chromosome
      top.marker <- mOnChr[which.max(marker.LRS[phe.name, mOnChr])]    # Marker showing the highest LRS on this chromosome
      mAboveT <- marker.LRS[phe.name, mOnChr] > LRS.threshold
      top.start <- names(which(mAboveT)[1])                            # Marker where the regions > LRS.threshold starts
      top.stop <- names(which(mAboveT)[length(which(mAboveT))])        # Marker where the regions > LRS.threshold stops

      strains <- rownames(phe.M.age.adj.SA)
      phenotype <- phe.M.age.adj.SA[strains, phe.name]
      marker.geno <- as.factor(genotypes[strains, top.marker])
      
      # Create a row in the output table for this QTL
      output.row <- c(paste0("QTL", qtl.idx), phe.name,
                      chrs[i], map[top.start, "Mb"], map[top.marker, "Mb"], map[top.stop, "Mb"],
                      top.marker, top.value,
                      mean(as.numeric(phenotype[which(marker.geno == "B")]), na.rm=TRUE),
                      mean(as.numeric(phenotype[which(marker.geno == "D")]), na.rm=TRUE),
                      sum(as.numeric(mAboveT)))
      output.tablular <- rbind(output.tablular, output.row)
      
      qtl.idx <- qtl.idx + 1
    }
  }
}
colnames(output.tablular) <- c("QTL ID", "Phenotype", "Chr", "Start", "Top", "Stop", "top marker", "max LRS", "Mean B", "Mean D", "M in region")
write.table(output.tablular, "output.tablular.femalesonly.SA.txt", sep = "\t", quote = FALSE, row.names=FALSE)





m <- "rs4226829"; phe <- "Tibia Ct.Th (mm) F"
strains <- rownames(phe.M.age.adj.SA)
phenotype <- phe.M.age.adj.SA[strains, phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]), notch = TRUE, xlab=m, ylab=phe, main="Raw data")





phenotype <- phe.M.age.adj.SA[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]), notch = TRUE, xlab=m, ylab=phe, main="Raw data")






phenotype <- "Femur Ct.pMOI (mm^4) F"
res <- NULL
for(marker in colnames(genotypes)) {
  res <- c(res, anova(lm(round(phe.M.age.adj.SA[,phenotype],3) ~ genotypes[,marker]))[[5]][1])
}
plot(-log10(res)*toLRS, t ='h')

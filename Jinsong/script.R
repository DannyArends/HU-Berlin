# Script by Danny Arends, for analysis of BXD bone density data by glm

library("lme4")
library("qtl")
library("ctl")

setwd("D:/Collegues/Jinsong")

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

# First, adjust for the major Sex effects we observe by using the residuals + mean
phe.sex.adj <- NULL
for(phe.name in phe.names){
  phe.adj <- get.residuals(get.eff(phenotypes, phe.name, "Sex")) +  mean(as.numeric(phenotypes[,phe.name]), na.rm = TRUE)
  phe.sex.adj <- cbind(phe.sex.adj, phe.adj)
}
colnames(phe.sex.adj) <- phe.names
phe.sex.adj <- cbind(phenotypes[,1:5], phe.sex.adj)
phe.sex.adj[1:5, 1:8]

# Then, adjust for the major logAge effects we observe by using the residuals + mean, on the sex adjusted data
phe.sex.age.adj <- NULL
for(phe.name in phe.names){
  phe.adj <- get.residuals(get.eff.num(phe.sex.adj, phe.name, "LogAGE")) +  mean(as.numeric(phe.sex.adj[,phe.name]), na.rm = TRUE)
  phe.sex.age.adj <- cbind(phe.sex.age.adj, phe.adj)
}
colnames(phe.sex.age.adj) <- phe.names
phe.sex.age.adj <- cbind(phenotypes[,1:5], phe.sex.age.adj)
phe.sex.age.adj[1:5, 1:8]

pheNoG <- which(!(phe.sex.age.adj[,"Strain"] %in% colnames(genotypes)))
genNoP <- which(!(colnames(genotypes) %in% phe.sex.age.adj[,"Strain"]))

# Keep the individuals that match between genotypes and phenotypes
phe.sex.age.adj <- phe.sex.age.adj[-pheNoG, ]
phe.sex.age.adj[1:5,1:5]

genotypes <- genotypes[, - genNoP]
genotypes <- t(genotypes[,phe.sex.age.adj[,"Strain"]])
rownames(genotypes) <- phe.sex.age.adj[,"Strain"]
genotypes[1:5,1:5]

dim(phe.sex.age.adj)
dim(genotypes)
toLRS <- 2*log(10)        # Conversion constant LOD -> LRS

# Resulting objects, LRS scoroes for the Sex:Genotype interaction, and the Genotype itself
marker.LRS      <- matrix(NA, length(phe.names), ncol(genotypes), dimnames=list(phe.names, colnames(genotypes)))
interaction.LRS <- matrix(NA, length(phe.names), ncol(genotypes), dimnames=list(phe.names, colnames(genotypes)))

for(phe.name in phe.names) {
  for(x in 1:ncol(genotypes)) {
    marker <- as.numeric(genotypes[,x] == "D") + 1                  # Numerical code for the marker, B = 1, D = 2
    phenotype <- as.numeric(phe.sex.age.adj[, phe.name])            # Phenotype as numeric values
    notNA <- which(!is.na(marker) & !is.na(phenotype))              # We can only use non-mising values
  
    marker.data <- marker[notNA];                                   #cat(phe.name, " @t ", colnames(genotypes)[x], ", m: ", length(marker.data), ", ")
    phenotype.data <- phenotype[notNA];                             #cat("p: ", length(phenotype.data), ", ")
    group.factor <- as.factor(phe.sex.age.adj[notNA, "Strain"]);    #cat("g: ", length(group.factor), ", ")
    sex.factor <- as.factor(phe.sex.age.adj[notNA, "Sex"]);         #cat("s: ", length(sex.factor), ", ")
    
    M.Full   <- lmer(phenotype.data ~ (1 | group.factor) + marker.data + marker.data:sex.factor, REML = FALSE)  # Full model
    M.Marker <- lmer(phenotype.data ~ (1 | group.factor) + marker.data, REML = FALSE)                           # Marker only model
    M.Null   <- lmer(phenotype.data ~ (1 | group.factor), REML = FALSE)                                         # Null model

    marker.LRS[phe.name, x] <- round(toLRS * -log10(anova(M.Marker, M.Null)$"Pr(>Chisq)"[2]), 2)                # Anova pvalue -> LOD -> toLRS for marker
    interaction.LRS[phe.name, x] <- round(toLRS * -log10(anova(M.Full, M.Marker)$"Pr(>Chisq)"[2]), 2)           # Anova pvalue -> LOD -> toLRS for interaction

    if(marker.LRS[phe.name, x] >= 12) {
      cat(phe.name, " @t ", colnames(genotypes)[x], ": ", marker.LRS[phe.name, x], ", ", interaction.LRS[phe.name, x], "\n")
    }
  }
}

# Write the results as marices to HDD
write.table(marker.LRS, "marker.LRS.txt", sep = "\t", quote = FALSE)
write.table(interaction.LRS, "interaction.LRS.txt", sep = "\t", quote = FALSE)

# Create plots of the different phenotypes
for(phe.name in phe.names) {
  plot(y = c(0, 25), x = c(0, max(map[,"cMplot"])), t = 'n', ylab = "LRS", xlab = "cM", main = phe.name)
  for(i in 1:length(chrs)){
    mOnChr <- rownames(map[which(map[,"Chr"] == chrs[i]),])
    chrID <- ((i%%2) + 1)
    points(marker.LRS[phe.name, mOnChr], x=map[mOnChr,"cMplot"], t = 'h', col=c("dodgerblue", "dodgerblue4")[chrID], lwd=1)
    points(interaction.LRS[phe.name, mOnChr], x=map[mOnChr,"cMplot"], t = 'h', col=c("burlywood4", "burlywood")[chrID], lwd=1)
  }
  abline(h = 15, col="green", lty = 2)
  abline(h = 12, col="gold", lty = 2)
}




get.eff.lme(phe.sex.age.adj, phe.name, "Strain")

lvls <- seq(0,1,length.out = nlevels(as.factor(phenotypes[,"Sex"])))
cols <- rgb(lvls, rep(1,length(lvls)), 1.0 - lvls, 0.5)
phe.adj.sex <- get.residuals(get.eff(phenotypes, phe.name, "Sex")) +  mean(as.numeric(phenotypes[,phe.name]), na.rm = TRUE)

op <- par(mfrow=c(2, 1))

boxplot(as.numeric(phenotypes[,phe.name]) ~ as.factor(phenotypes[,"Sex"]), col = cols, notch = TRUE)
boxplot(phe.adj.sex ~ as.factor(phenotypes[,"Sex"]), add=TRUE, col = cols, notch = TRUE)

phe.adj.logAge <- get.residuals(get.eff.num(phenotypes, phe.name, "LogAGE")) +  mean(as.numeric(phenotypes[,phe.name]), na.rm = TRUE)

plot(as.numeric(phenotypes[,phe.name]) ~ as.numeric(phenotypes[,"LogAGE"]), col = rgb(1.0, 0.0, 0.0, 0.7), pch = 18, ylab = phe.name, main = phe.name)
points(phe.adj.logAge ~ as.numeric(phenotypes[,"LogAGE"]), col = rgb(0.5, 1.0, 0.2, 0.7), pch = 18)
legend("topright", c("Before Age correction", "After Age correction") ,col =c(rgb(1.0, 0.0, 0.0, 0.7), rgb(0.5, 1.0, 0.2, 0.7)), pch = 18, bg  = "white")




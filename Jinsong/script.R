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

# Write the results as matrices to HDD
write.table(marker.LRS, "marker.LRS.txt", sep = "\t", quote = FALSE)
write.table(interaction.LRS, "interaction.LRS.txt", sep = "\t", quote = FALSE)

marker.LRS <- read.table("marker.LRS.txt", sep = "\t",check.names=FALSE)
interaction.LRS <- read.table("interaction.LRS.txt", sep = "\t",check.names=FALSE)

# Create plots of the different phenotypes
for(phe.name in phe.names) {
  png(paste0(gsub("%", "Pct", gsub("/", ".", gsub(" F", "", phe.name))), ".png"), width = 1200, height =  600)
  plot(y = c(0, 25), x = c(0, max(map[,"cMplot"])), t = 'n', ylab = "LRS", xlab = "cM", main = gsub(" F", "", phe.name))
  for(i in 1:length(chrs)){
    mOnChr <- rownames(map[which(map[,"Chr"] == chrs[i]),])
    chrID <- ((i%%2) + 1)
    points(marker.LRS[phe.name, mOnChr], x=map[mOnChr,"cMplot"], t = 'h', col=c("dodgerblue", "dodgerblue4")[chrID], lwd=1)
    points(interaction.LRS[phe.name, mOnChr], x=map[mOnChr,"cMplot"], t = 'h', col=c("burlywood4", "burlywood")[chrID], lwd=1)
  }
  abline(h = 15, col="green", lty = 2)
  abline(h = 12, col="gold", lty = 2)
  dev.off()
}

LRS.cutoff <- 15.00
LRS.threshold <- LRS.cutoff - (1.5 * toLRS)

# Create an output table for QTL results
output.tablular <- NULL
qtl.idx <- 1
for(phe.name in phe.names) {
  for(i in 1:length(chrs)) {
    mOnChr <- rownames(map[which(map[,"Chr"] == chrs[i]),])         # Markers on chromosome
    mAboveC <- marker.LRS[phe.name, mOnChr] >= LRS.cutoff            # Markers above 13.5 LRS
    if(any(mAboveC)){
      top.value <- max(marker.LRS[phe.name, mOnChr], na.rm = TRUE)     # Top LRS score on this chromosome
      top.marker <- mOnChr[which.max(marker.LRS[phe.name, mOnChr])]    # Marker showing the highest LRS on this chromosome
      mAboveT <- marker.LRS[phe.name, mOnChr] >= LRS.threshold
      top.start <- mOnChr[which(mAboveT)][1]                            # Marker where the regions > LRS.threshold starts
      top.stop <- mOnChr[which(mAboveT)[length(which(mAboveT))]]        # Marker where the regions > LRS.threshold stops

      # For the plot we need the phenotype data, and the genotype data at the top.marker
      strains <- phe.sex.age.adj[which(phe.sex.age.adj[,"Strain"] %in% rownames(genotypes)),"Strain"]
      phenotype <- phe.sex.age.adj[which(phe.sex.age.adj[,"Strain"] %in% rownames(genotypes)), phe.name]
      marker.geno <- as.factor(genotypes[strains, top.marker])
      
      # Create an effect plot on the corrected phenotype data for each QTL
      png(paste0("QTL", qtl.idx, ".png"), width=400, height=400)
        chr.pos <- paste0(chrs[i],":", map[top.start, "Mb"], "-", map[top.stop, "Mb"])
        plot(as.numeric(phenotype) ~ marker.geno, notch = TRUE, xlab="", ylab=phe.name, main = paste0(top.marker," (", chr.pos,")"))
      dev.off()

      # Create a row in the output table for this QTL
      output.row <- c(paste0("QTL", qtl.idx), phe.name,
                      chrs[i], map[top.start, "Mb"], map[top.marker, "Mb"], map[top.stop, "Mb"],
                      top.marker, top.value, interaction.LRS[phe.name, top.marker],
                      mean(as.numeric(phenotype[which(marker.geno == "B")]), na.rm=TRUE),
                      mean(as.numeric(phenotype[which(marker.geno == "D")]), na.rm=TRUE),
                      sum(as.numeric(mAboveT)))
      output.tablular <- rbind(output.tablular, output.row)
      qtl.idx <- qtl.idx + 1
    }
  }
}
colnames(output.tablular) <- c("QTL ID", "Phenotype", "Chr", "Start", "Top", "Stop", "top marker", "max LRS", "G:SEX LRS", "Mean B", "Mean D", "M in region")

# Write the output table with QTL results to HDD
write.table(output.tablular, "output.tablular.txt", sep = "\t", quote = FALSE, row.names=FALSE)

output.tablular <- read.table("output.tablular.txt", sep = "\t", header=TRUE)

# Double check for G:Sex effects
for(phe.name in phe.names) {
  if(any(interaction.LRS[phe.name, ] > LRS.cutoff)){
    cat(phe.name, "might have a strong G:Sex effect", "\n")
  }
}

# Show 4 QTLs in the RAW data (LRS above 15, should be visible in raw data more or less)
strains <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)),"Strain"]
sex     <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), "Sex"]

op <- par(mfrow=c(1, 2))

m <- "rs3659436"; phe <- "Tibia Trab.DA (ratio) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Raw data")

m <- "rs3659436"; phe <- "Tibia Trab.DA (ratio) F"
phenotype <- phe.sex.age.adj[which(phe.sex.age.adj[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Adj data")



m <- "rs3653769"; phe <- "Tibia Ct.BV (mm^3) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ sex + as.factor(genotypes[strains, m]), notch = TRUE, xlab=m, ylab=phe, main="Raw data")

m <- "rs13478223"; phe <- "Femur Ct.Apparent.BMD (mgHA/ cm^3) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ sex + as.factor(genotypes[strains, m]), notch = TRUE, xlab=m, ylab=phe, main="Raw data")

m <- "mCV25103990"; phe <- "Tibia Ct. Porosity (%) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ sex + as.factor(genotypes[strains, m]), notch = TRUE, xlab=m, ylab=phe, main="Raw data")

# We observe that the QTLs we find here, have no sex specific component
# Directions of effect are the same between M.B and M.D versus F.B and F.D


### Comparison to Jinsong

m <- "rs3682996"; phe <- "Tibia Trab TV (mm^3) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Raw phenotype data (Fttf1a)", sub = "Why is the LRS for males 1 in GN?")

m <- "rs3687595"; phe <- "Femur Ct.Th (mm) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Raw phenotype data (Fttf1b)")


op <- par(mfrow=c(1,2))

m <- "rs3657682"; phe <- "Femur Ct.TV (mm^3) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Raw phenotype data (Fcvf12)")


m <- "rs3657682"; phe <- "Femur Ct.TV (mm^3) F"
phenotype <- phe.sex.adj[which(phe.sex.adj[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="phenotype data sex adj (Fcvf12)")


m <- "rs3682996"; phe <- "Femur Ct.TV (mm^3) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Raw phenotype data (Fcvf12)")






op <- par(mfrow=c(1, 2))

m <- "CEL-8_25820220"; phe <- "Femur Length (mm) F"
phenotype <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Raw data")

m <- "CEL-8_25820220"; phe <- "Femur Length (mm) F"
phenotype <- phe.sex.age.adj[which(phe.sex.age.adj[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Adj data")





# Show 4 QTLs in the RAW data (LRS above 15, should be visible in raw data more or less)
strains <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)),"Strain"]
sex     <- phenotypes[which(phenotypes[,"Strain"] %in% rownames(genotypes)), "Sex"]




m <- "rs4226829"; phe <- "Tibia Ct.Th (mm) F"
phenotype <- phe.sex.age.adj[which(phe.sex.age.adj[,"Strain"] %in% rownames(genotypes)), phe]
boxplot(as.numeric(phenotype) ~ as.factor(genotypes[strains, m]) + sex, notch = TRUE, xlab=m, ylab=phe, main="Adj data")




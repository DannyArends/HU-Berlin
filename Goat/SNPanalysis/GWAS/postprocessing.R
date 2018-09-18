setwd("D:/Edrive/Goat/DNA/SihamAnalysis/GWAS")
genotypes <- read.table("genotypes.txt", sep="\t", check.names=FALSE)
map <- read.table("geneticmap.txt", sep="\t", check.names=FALSE)
samples <- read.table("samples.txt", sep="\t", check.names=FALSE)

dim(genotypes)
genotypes[1:5,1:5]
dim(map)
map[1:5,]
dim(samples)
samples[1:5,]

pheNames <- c("Averagemilk", "Weight", "Withersheight", "Rumpheight", "Bodylength", "Sternumheight", 
              "Bodydepth", "Bicoastaldiameter", "Earlength", "RumpWidth", "HeadWidth", "Rumplength", 
              "Headlength", "Heartgirth", "Cannonbone", "Muzzlediameter")

breeds <- c("Dese", "Ni", "Nu", "Tagg")
              
phenotypes <- apply(samples[, pheNames], 2, as.numeric)

# No milk, since they are not lactating yet
noMilk <- which(phenotypes[,"Averagemilk"] == 0)
phenotypes[noMilk,"Averagemilk"] <- NA

# Breeds of the different animals
breedV <- samples[,"Breed"]

phenotypes.adjusted <- phenotypes
rownames(phenotypes.adjusted) <- 1:nrow(phenotypes.adjusted)

# Correct for Age and Breed
adjustments <- NULL
for(pheno in pheNames){
  mymodel_A  <- lm(phenotypes[,pheno] ~ log(as.numeric(samples[,"Age"])))
  mymodel_B  <- lm(phenotypes[,pheno] ~ breedV)
  mymodel_AB <- lm(phenotypes[,pheno] ~ log(as.numeric(samples[,"Age"])) + breedV)
  #cat(pheno, anova(mymodel_AB)[[5]][1], anova(mymodel_AB)[[5]][2], anova(mymodel_AB)[[5]][3], "\n")
  # We can have 3 classes: Age affected, Breed affected and Age and Breed affected
  if(anova(mymodel_AB)[[5]][1] < 0.1 && anova(mymodel_AB)[[5]][2] < 0.1){
    mymodel <- mymodel_AB
    cat("Adjusting:", pheno, " Age and Breed\n")
    phenoresiduals <- mymodel$residuals
    phenotypes.adjusted[names(phenoresiduals), pheno] <- mean(phenotypes[,pheno],na.rm = TRUE) + mymodel$residuals
  }else if(anova(mymodel_AB)[[5]][1] < 0.1 && anova(mymodel_AB)[[5]][2] >= 0.1){
    mymodel <- mymodel_A
    cat("Adjusting:", pheno, " Age\n")
    phenoresiduals <- mymodel$residuals
    phenotypes.adjusted[names(phenoresiduals), pheno] <- mean(phenotypes[,pheno],na.rm = TRUE) + mymodel$residuals
  }else if(anova(mymodel_AB)[[5]][1] >= 0.1 && anova(mymodel_AB)[[5]][2] < 0.1){
    mymodel <- mymodel_B
    cat("Adjusting:", pheno, " Breed\n")
    phenoresiduals <- mymodel$residuals
    phenotypes.adjusted[names(phenoresiduals), pheno] <- mean(phenotypes[,pheno],na.rm = TRUE) + mymodel$residuals
  }else{
    # No adjustment
  }
  adjustments <- rbind(adjustments, c(paste0(round(mean(phenotypes[breedV == "Dese", pheno], na.rm = TRUE),1), " (",round(sd(phenotypes[breedV == "Dese", pheno], na.rm = TRUE),1), ")"), 
                                      paste0(round(mean(phenotypes[breedV == "Ni", pheno], na.rm = TRUE),1), " (",round(sd(phenotypes[breedV == "Ni", pheno], na.rm = TRUE),1), ")"),
                                      paste0(round(mean(phenotypes[breedV == "Nu", pheno], na.rm = TRUE),1), " (",round(sd(phenotypes[breedV == "Nu", pheno], na.rm = TRUE),1), ")"),
                                      paste0(round(mean(phenotypes[breedV == "Tagg", pheno], na.rm = TRUE),1), " (",round(sd(phenotypes[breedV == "Tagg", pheno], na.rm = TRUE),1), ")"), anova(mymodel_AB)[[5]][1:2]))
}
rownames(phenotypes.adjusted) <- rownames(phenotypes)
rownames(adjustments) <- colnames(phenotypes)
colnames(adjustments) <- c("Desert", "Nilotic", "Nubian", "Taggar", "P(age)", "P(breed)")

write.table(adjustments, "adjustments.txt", sep="\t", quote=FALSE)

# Find outliers, and put them to missing values
for(pheno in pheNames){
  above <- which(phenotypes.adjusted[,pheno] > mean(phenotypes.adjusted[,pheno],na.rm = TRUE) + 3 * sd(phenotypes.adjusted[,pheno], na.rm=TRUE))
  below <- which(phenotypes.adjusted[,pheno] < mean(phenotypes.adjusted[,pheno],na.rm = TRUE) - 3 * sd(phenotypes.adjusted[,pheno], na.rm=TRUE))
  if(length(above) > 0 || length(below) > 0) {
    cat("Outlier in", pheno, length(above), length(below), "\n")
    phenotypes.adjusted[c(above,below), pheno] <- NA
  }
}

phenotypes.adjusted <- phenotypes.adjusted[,-c(1:2)]

setwd("D:/Edrive/Goat/DNA/SihamAnalysis/GWAS")
lod.add <- read.table("lod.add.txt", sep="\t")
pvalues.add <- read.table("pvalues.add.txt", sep="\t")
lod.dom <- read.table("lod.dom.txt", sep="\t")
pvalues.dom <- read.table("pvalues.dom.txt", sep="\t")
map <- read.table("geneticmap.txt", sep="\t", check.names=FALSE)


pheNames <- c("Averagemilk", "Weight", "Withersheight", "Rumpheight", "Bodylength", "Sternumheight", 
              "Bodydepth", "Bicoastaldiameter", "Earlength", "RumpWidth", "HeadWidth", "Rumplength", 
              "Headlength", "Heartgirth", "Cannonbone", "Muzzlediameter")

# Get the significant and suggestive amounts
lod.add.adj <- lod.add
pvalues.add.adj <- pvalues.add

lod.dom.adj <- lod.dom
pvalues.dom.adj <- pvalues.dom


genotypes["snp2253-scaffold1069-1041798",]

phenotype <- "Bicoastaldiameter"
sampleSize <- table(as.character(unlist(genotypes[maxMarker,])))
maxMarker <- rownames(lod.add.adj)[which.max(lod.add.adj[,phenotype])]
m <- as.numeric(as.factor(as.character(unlist(genotypes[maxMarker,]))))
model <- lm(as.numeric(phenotypes.adjusted[,phenotype]) ~ m)

n1 <-  phenotypes.adjusted[which(as.character(unlist(genotypes[maxMarker,])) == "AA"),phenotype]
n2 <-  phenotypes.adjusted[which(as.character(unlist(genotypes[maxMarker,])) == "AG"),phenotype]
n3 <-  phenotypes.adjusted[which(as.character(unlist(genotypes[maxMarker,])) == "GG"),phenotype]
p1 <- t.test(n1, n2)$p.value# * (3 * 24027)
p2 <- t.test(n1, n3)$p.value# * (3 * 24027)
p3 <- t.test(n2, n3)$p.value# * (3 * 24027)

par(cex = 1.2)
plot(c(0.5, 3.5), c(5, 45), t ='n', main=paste0(phenotype," @ ",maxMarker), xlab="", ylab=paste0(phenotype, " (cm)"), xaxt='n')
par(cex = 1)
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), varwidth = TRUE, add = TRUE,xaxt='n', yaxt='n')
par(cex = 1.2)
axis(1, at = c(1,2,3), paste0(names(sampleSize), " (N = ", sampleSize, ")"))
lines(c(1, 2) , y = c(33,33))
lines(c(1, 1) , y = c(32,33))
lines(c(2, 2) , y = c(32,33))
text(1.5, 34, "**")
lines(c(2, 3) , y = c(37,37))
lines(c(2, 2) , y = c(36,37))
lines(c(3, 3) , y = c(36,37))
text(2.5, 38, "**")
lines(c(1, 3) , y = c(41,41))
lines(c(1, 1) , y = c(40,41))
lines(c(3, 3) , y = c(40,41))
text(2, 42, "***")

phenotype <- "Bodylength"
sampleSize <- table(as.character(unlist(genotypes[maxMarker,])))
maxMarker <- rownames(lod.add.adj)[which.max(lod.add.adj[,phenotype])]

n1 <-  phenotypes.adjusted[which(as.character(unlist(genotypes[maxMarker,])) == "CC"),phenotype]
n2 <-  phenotypes.adjusted[which(as.character(unlist(genotypes[maxMarker,])) == "TC"),phenotype]
n3 <-  phenotypes.adjusted[which(as.character(unlist(genotypes[maxMarker,])) == "TT"),phenotype]
p1 <- t.test(n1, n2)$p.value# * (3 * 24027)
p2 <- t.test(n1, n3)$p.value# * (3 * 24027)
p3 <- t.test(n2, n3)$p.value# * (3 * 24027)


par(cex = 1.2)
plot(c(0.5, 3.5), c(45, 80), t ='n', main=paste0(phenotype," @ ",maxMarker), xlab="", ylab=paste0(phenotype, " (cm)"), xaxt='n')
par(cex = 1)
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), varwidth = TRUE, add = TRUE,xaxt='n', yaxt='n')
par(cex = 1.2)
axis(1, at = c(1,2,3), paste0(names(sampleSize), " (N = ", sampleSize, ")"))
lines(c(1, 2) , y = c(75,75))
lines(c(1, 1) , y = c(74,75))
lines(c(2, 2) , y = c(74,75))
text(1.5, 76, "**")
lines(c(2, 3) , y = c(72,72))
lines(c(2, 2) , y = c(71,72))
lines(c(3, 3) , y = c(71,72))
text(2.5, 73, "**")
lines(c(1, 3) , y = c(79,79))
lines(c(1, 1) , y = c(78,79))
lines(c(3, 3) , y = c(78,79))
text(2, 80, "***")


abline(model$coefficients["(Intercept)"], model$coefficients["m"])

as.numeric(phenotypes.adjusted[,phenotype]) 

op <- par(mfrow=c(1,5))
phenotype <- "HeadWidth"
maxMarker <- rownames(lod.add.adj)[which.max(lod.add.adj[,phenotype])]
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), varwidth = TRUE, main=paste0(phenotype," @ ",maxMarker), xlab="", ylab=paste0(phenotype, " (cm)"))

phenotype <- "Rumplength"
maxMarker <- rownames(lod.add.adj)[which.max(lod.add.adj[,phenotype])]
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), varwidth = TRUE, main=paste0(phenotype," @ ",maxMarker), xlab="", ylab=paste0(phenotype, " (cm)"))


phenotype <- "Rumplength"
maxMarker <- "snp33453-scaffold392-4238920"
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), varwidth = TRUE, main=paste0(phenotype," @ ",maxMarker), xlab="", ylab=paste0(phenotype, " (cm)"))

phenotype <- "Cannonbone"
maxMarker <- rownames(lod.dom.adj)[which.max(lod.dom.adj[,phenotype])]
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), varwidth = TRUE, main=paste0(phenotype," @ ",maxMarker), xlab="", ylab=paste0(phenotype, " (cm)"))


phenotype <- "Withersheight"
maxMarker <- rownames(lod.dom.adj)[which.max(lod.dom.adj[,phenotype])]
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), varwidth = TRUE, main=paste0(phenotype," @ ",maxMarker), xlab="", ylab=paste0(phenotype, " (cm)"))


lod.add.adj["snp44764-scaffold609-1510371", "Rumplength"]
lod.dom.adj["snp44764-scaffold609-1510371", "Rumplength"]
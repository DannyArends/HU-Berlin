#
# Goat GWAS
#

createAddDom <- function(geno.char){
  alleles <- sort(unique(unlist(strsplit(geno.char, ""))))
  geno.add <- rep(NA, length(geno.char))
  geno.dom <- rep(NA, length(geno.char))
  AA <- paste0(alleles[1],alleles[1])
  AB <- paste0(alleles[1],alleles[2])
  BA <- paste0(alleles[2],alleles[1])
  BB <- paste0(alleles[2],alleles[2])
  geno.add[geno.char == AA] <- 0
  geno.add[geno.char == AB] <- 1
  geno.add[geno.char == BA] <- 1
  geno.add[geno.char == BB] <- 2

  geno.dom[geno.char == AA] <- 0
  geno.dom[geno.char == AB] <- 1
  geno.dom[geno.char == BA] <- 1
  geno.dom[geno.char == BB] <- 0
  return(list(geno.add = geno.add, geno.dom = geno.dom))
}

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

pvalues.add <- matrix(NA, nrow(genotypes), length(pheNames), dimnames=list(rownames(genotypes), pheNames))
pvalues.dom <- matrix(NA, nrow(genotypes), length(pheNames), dimnames=list(rownames(genotypes), pheNames))
cnt <- 1
for(marker in rownames(genotypes)){
  geno.char <- as.character(unlist(genotypes[marker,]))
  recoded <- createAddDom(geno.char)
  for(pheno in pheNames){
    res <- anova(lm(phenotypes.adjusted[,pheno] ~ recoded$geno.add + recoded$geno.dom))
    pvals <- na.omit(res$"Pr(>F)")
    if(length(pvals) == 2){
      pvalues.add[marker, pheno] <- pvals[1]
      pvalues.dom[marker, pheno] <- pvals[2]
    }
  }
  cat("Done", cnt, "/", nrow(genotypes), "\n")
  cnt <- cnt + 1
}

setwd("D:/Edrive/Goat/DNA/SihamAnalysis/GWAS")
write.table(pvalues.add, "pvalues.add.txt", sep="\t", quote=FALSE)
write.table(pvalues.dom, "pvalues.dom.txt", sep="\t", quote=FALSE)

lod.add <- round(-log10(pvalues.add),2)
lod.dom <- round(-log10(pvalues.dom),2)

write.table(lod.add, "lod.add.txt", sep="\t", quote=FALSE)
write.table(lod.dom, "lod.dom.txt", sep="\t", quote=FALSE)

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

signsug <- NULL
for(pheno in pheNames){
  observed.p.add <- pvalues.add[, pheno]
  expected.p.add <- (rank(observed.p.add, ties.method="first")+0.5) / (length(observed.p.add) + 1)
  expected.l.add <- -log10(sort(expected.p.add, na.last = TRUE))
  observed.l.add <- -log10(sort(observed.p.add, na.last = TRUE))

  observed.p.dom <- pvalues.dom[, pheno]
  expected.p.dom <- (rank(observed.p.dom, ties.method="first")+0.5) / (length(observed.p.dom) + 1)
  expected.l.dom <- -log10(sort(expected.p.dom, na.last = TRUE))
  observed.l.dom <- -log10(sort(observed.p.dom, na.last = TRUE))

  # Estimate the lambda value
  lambda.add <- round(median(qchisq(1.0 - observed.p.add, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)
  lambda.dom <- round(median(qchisq(1.0 - observed.p.dom, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)
  
  zscores.add <- qnorm(1 - (observed.p.add/2)) 
  zscores.dom <- qnorm(1 - (observed.p.dom/2)) 
  
  # Adjustment of Pvalues based on the lambda
  if(diff(c(1,lambda.add)) > 0.05) {
    cat("Adjusting lambda ")
    observed.p.add.adj = pchisq((zscores.add^2)/lambda.add, df = 1, lower = FALSE)
    expected.p.add.adj <- (rank(observed.p.add.adj, ties.method="first")+0.5) / (length(observed.p.add.adj) + 1)
    expected.l.add.adj <- -log10(sort(expected.p.add.adj, na.last = TRUE))
    observed.l.add.adj <- -log10(sort(observed.p.add.adj, na.last = TRUE))
    pvalues.add.adj[,pheno] <- observed.p.add.adj
    lod.add.adj[,pheno] <- -log10(observed.p.add.adj)
  }else{
    observed.p.add.adj = observed.p.add
    expected.p.add.adj = expected.p.add
    observed.l.add.adj = observed.l.add
    expected.l.add.adj = expected.l.add
  }
  
  if(diff(c(1,lambda.dom)) > 0.05) {
    cat("Adjusting lambda ")
    observed.p.dom.adj = pchisq((zscores.dom^2)/lambda.dom, df = 1, lower = FALSE)
    expected.p.dom.adj <- (rank(observed.p.dom.adj, ties.method="first")+0.5) / (length(observed.p.dom.adj) + 1)
    expected.l.dom.adj <- -log10(sort(expected.p.dom.adj, na.last = TRUE))
    observed.l.dom.adj <- -log10(sort(observed.p.dom.adj, na.last = TRUE))
    pvalues.dom.adj[,pheno] <- observed.p.dom.adj
    lod.dom.adj[,pheno] <- -log10(observed.p.dom.adj)
  }else{
    observed.p.dom.adj = observed.p.dom
    expected.p.dom.adj = expected.p.dom
    observed.l.dom.adj = observed.l.dom
    expected.l.dom.adj = expected.l.dom
  }

  # Adjusted lambda (should be 1)
  lambda.add.adj <- round(median(qchisq(1.0 - observed.p.add.adj, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)
  lambda.dom.adj <- round(median(qchisq(1.0 - observed.p.dom.adj, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)

  # Create the plot
  png(paste0("QQplots/QQplot",pheno,".png"))
    plot(c(0,8), c(0,8), main = paste0("QQ plot: ", pheno), t ='n')
    points(expected.l.add, observed.l.add, pch=19,cex=0.8)
    points(expected.l.add.adj, observed.l.add.adj, col="green", pch=19,cex=0.7)
  
    points(expected.l.dom, observed.l.dom, col="gray", pch=17,cex=0.8)
    points(expected.l.dom.adj, observed.l.dom.adj, col="blue", pch=17,cex=0.7)
  
    abline(0,1)
    legend("topleft", c("Additive", "Additive (Adj)", "DomDev", "DomDev (Adj)"), col=c("black","green", "gray", "blue"), pch=c(19,19,17,17))
  dev.off()
  cat(pheno, lambda.add, lambda.dom, lambda.add.adj, lambda.dom.adj, "\n")

  # Find significant/suggestive effects
  signlvl <- 5.6
  
  sign.add <- length(which(observed.l.add > signlvl));sug.add <- length(which(observed.l.add > 5))
  sign.dom <- length(which(observed.l.dom > signlvl));sug.dom <- length(which(observed.l.dom > 5))

  sign.add.adj <- length(which(observed.l.add.adj > signlvl));sug.add.adj <- length(which(observed.l.add.adj > 5))
  sign.dom.adj <- length(which(observed.l.dom.adj > signlvl));sug.dom.adj <- length(which(observed.l.dom.adj > 5))

  signsug <- rbind(signsug, c(sign.add, sug.add, sign.dom, sug.dom, sign.add.adj, sug.add.adj, sign.dom.adj, sug.dom.adj))
}

colnames(signsug) <- c("sign.add", "sug.add", "sign.dom", "sug.dom", "sign.add.adj", "sug.add.adj", "sign.dom.adj", "sug.dom.adj")
rownames(signsug) <- pheNames

chrs <- c(1:29,"X")
gap <- 20000000
map.sorted <- NULL
chr.lengths <- c()
chr.starts <- c(0)
i <- 1
for(chr in chrs){
  onChr <- which(map[,"Chr"] == chr)
  map.sorted <- rbind(map.sorted, map[onChr,])
  chr.lengths <- c(chr.lengths, max(map[onChr, "Position"]))
  chr.starts <- c(chr.starts, chr.starts[i] + max(map[onChr, "Position"]) + gap)
  i <- i + 1
}

for(pheno in pheNames){
  png(paste0("GWASplots/GWAS",pheno,"_Add.png"), width=1024)
    plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", pheno, ", Additive Effect"))
    for(chr in chrs){
      onChr <- rownames(map.sorted[map.sorted[,"Chr"] == chr,])
      points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = -log10(pvalues.add.adj)[onChr, pheno],t ='h', col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
      if(chr != "X") text(x=chr.starts[chr] + chr.lengths[chr] / 2, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
      if(chr == "X") text(x=chr.starts[chr] + 2*gap, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
      i <- i + 1
    }
    abline(h=-log10(0.05/24027), col="green",lty=3)
    abline(h=5, col="orange",lty=3)
  dev.off()
  png(paste0("GWASplots/GWAS",pheno,"_DomDev.png"), width=1024)
    plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", pheno, ", Dominance Deviation"))
    for(chr in chrs){
      onChr <- rownames(map.sorted[map.sorted[,"Chr"] == chr,])
      points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = -log10(pvalues.dom.adj)[onChr, pheno],t ='h', col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
      if(chr != "X") text(x=chr.starts[chr] + chr.lengths[chr] / 2, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
      if(chr == "X") text(x=chr.starts[chr] + 2*gap, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
      i <- i + 1
    }
    abline(h=-log10(0.05/24027), col="green",lty=3)
    abline(h=5, col="orange",lty=3)
  dev.off()
}


names(chr.starts) <- chrs
names(chr.lengths) <- chrs
phenotype <- "Bicoastaldiameter"

layout(matrix(c(1,1,2,3),2,2,byrow=T))
op <- par(mar = c(4,4,2,1))
i <- 1
maxMarker <- rownames(lod.add.adj)[which.max(lod.add.adj[,phenotype])]

map[maxMarker,]
lod.add.adj[which.max(lod.add.adj[,phenotype]),phenotype]

plot(x = c(-gap, tail(chr.starts,1)), y = c(0,8), t = 'n', xlab="", ylab="-log10(P)",xaxt='n', xaxs="i", yaxs="i",las=2,main=paste0("Manhattan plot - ", phenotype, ", Additive Effect"))
for(chr in chrs){
  onChr <- rownames(map.sorted[map.sorted[,"Chr"] == chr,])
  points(x=chr.starts[chr] + map.sorted[onChr,"Position"], y = -log10(pvalues.add.adj)[onChr, phenotype],t ='h', col=c("black", "cornflowerblue")[(i %% 2 == 1) + 1])
  if(chr != "X") text(x=chr.starts[chr] + chr.lengths[chr] / 2, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
  if(chr == "X") text(x=chr.starts[chr] + 2*gap, y = 0.25, chr, col=c("white", "black")[(i %% 2 == 1) + 1],cex=0.8)
  i <- i + 1
}
abline(h=-log10(0.05/24027), col="green",lty=3)
abline(h=5, col="orange",lty=3)

observed.p.ori <- pvalues.add[, phenotype]
expected.p.ori <- (rank(observed.p.ori, ties.method="first")+0.5) / (length(observed.p.ori) + 1)

expected.l.ori <- -log10(sort(expected.p.ori, na.last = TRUE))
observed.l.ori <- -log10(sort(observed.p.ori, na.last = TRUE))

# Estimate the lambda value
lambda.ori <- round(median(qchisq(1.0 - observed.p.ori, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)


observed.p <- pvalues.add.adj[, phenotype]
expected.p <- (rank(observed.p, ties.method="first")+0.5) / (length(observed.p) + 1)

expected.l <- -log10(sort(expected.p, na.last = TRUE))
observed.l <- -log10(sort(observed.p, na.last = TRUE))

# Estimate the adjusted lambda value
lambda <- round(median(qchisq(1.0 - observed.p, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)
cat(phe, lambda, "\n")
#}

plot(c(0,4.5),c(0,9), main=paste0("QQ plot of ", phenotype), xlab="Expected", ylab="Observed", t = 'n')
points(expected.l, observed.l, pch=19, cex=0.8)
points(expected.l.ori, observed.l.ori, pch=19, cex=0.8, col="cornflowerblue")
legend("topleft", c(paste0("Original λ = ", lambda.ori), "Adjusted λ = 1.00"), pch=19, cex=0.8, col=c("cornflowerblue","black"))


abline(0, 1)
boxplot(as.numeric(phenotypes.adjusted[,phenotype]) ~ as.character(unlist(genotypes[maxMarker,])), main=paste0(phenotype," @ ",maxMarker), xlab="Genotype", ylab=paste0(phenotype, " (cm)"),varwidth = TRUE)



marker <- which.max(lod.add[,phenotype])
geno.char <- as.character(unlist(genotypes[marker,]))
recoded <- createAddDom(geno.char)
boxplot(phenotypes.adjusted[,phenotype] ~ recoded$geno.add, main=phenotype)

### Protein coding genes in the neaghbourhood

map <- read.table("../FilteredLocationLWLT01.txt", sep="\t", check.names=FALSE, header=TRUE, row.names=1)
map2 <- read.table("../FilteredLocationChir1.0.txt", sep="\t", check.names=FALSE, header=TRUE, row.names=1)
map3 <- read.table("../FilteredLocationChir2.0.txt", sep="\t", check.names=FALSE, header=TRUE, row.names=1)
map <- cbind(map, Position = (map[,"Start"] + map[,"Stop"]) / 2)
proteins <- read.csv("../../annotation/LWLT01/GeneList_01Jun2017.txt",sep="\t")

#proteins[which(as.character(proteins[,1]) == 2 & 
#                             proteins[,3] >  76820044 - 1000000 & 
#                             proteins[,3] <  76820044 + 1000000),]

                             
lod.dom.adj <- round(-log10(pvalues.dom.adj),2)
lod.add.adj <- round(-log10(pvalues.add.adj),2)

for(x in 1: ncol(lod.dom.adj)){
  i.dom <- which(lod.dom.adj[, x] > 5.6)
  i.add <- which(lod.add.adj[, x] > 5.6)
  if(length(i.dom) > 0){
    cat(x, "sign dom for", colnames(lod.dom.adj)[x], rownames(lod.dom.adj)[i.dom], "\n")
    chr <- map[i.dom,"chrN"]
    pos <- as.numeric(map[i.dom,"Position"])
    cat(chr, pos, "\n")
    write.table(proteins[which(as.character(proteins[,"chromosome"]) == chr & 
                             proteins[,"end_position_on_the_genomic_accession"] >  pos - 2000000 & 
                             proteins[,"start_position_on_the_genomic_accession"] <  pos + 2000000),], file = paste0("proteins_sign_dom_near_",colnames(lod.dom.adj)[x],".txt"), sep="\t", quote=FALSE)
  }
  if(length(i.add) > 0){
    cat(x, "sign dom for", colnames(lod.dom.adj)[x], rownames(lod.dom.adj)[i.add], "\n")
    probe <- rownames(lod.dom.adj)[i.add]
    chr <- map[probe,"chrN"]
    pos <- as.numeric(map[probe,"Position"])
    cat(chr, pos, "\n")
    write.table(proteins[which(as.character(proteins[,"chromosome"]) == chr & 
                             proteins[,"end_position_on_the_genomic_accession"] >  pos - 2000000 & 
                             proteins[,"start_position_on_the_genomic_accession"] <  pos + 2000000),], file = paste0("proteins_sign_add_near_",colnames(lod.add.adj)[x],".txt"), sep="\t", quote=FALSE)
  }
}


for(x in 1: ncol(lod.dom.adj)){
  ii.dom <- which(lod.dom.adj[, x] > 5)
  ii.add <- which(lod.add.adj[, x] > 5)
  if(length(ii.dom) > 0){
    for(i.dom in ii.dom){
    cat(x, "sign dom for", colnames(lod.dom.adj)[x], rownames(lod.dom.adj)[i.dom], "\n")
    probe <- rownames(lod.dom.adj)[i.dom]
    chr <- map[probe,"chrN"]
    pos <- as.numeric(map[probe,"Position"])
    cat(chr, pos, "\n")
    write.table(proteins[which(as.character(proteins[,"chromosome"]) == chr & 
                             proteins[,"end_position_on_the_genomic_accession"] >  pos - 2000000 & 
                             proteins[,"start_position_on_the_genomic_accession"] <  pos + 2000000),], file = paste0("gene_dom_near_",probe,"_",colnames(lod.dom.adj)[x],".txt"), sep="\t", quote=FALSE)
    }
  }
  if(length(ii.add) > 0){
    for(i.add in ii.add){
      cat(x, "sign add for", colnames(lod.dom.adj)[x], rownames(lod.dom.adj)[i.add], "\n")
      probe <- rownames(lod.dom.adj)[i.add]
      chr <- map[probe,"chrN"]
      pos <- as.numeric(map[probe,"Position"])
      cat(chr, pos, "\n")
      write.table(proteins[which(as.character(proteins[,"chromosome"]) == chr & 
                             proteins[,"end_position_on_the_genomic_accession"] >  pos - 2000000 & 
                             proteins[,"start_position_on_the_genomic_accession"] <  pos + 2000000),], file = paste0("gene_add_near_",probe,"_",colnames(lod.add.adj)[x],".txt"), sep="\t", quote=FALSE)
    }
  }
}

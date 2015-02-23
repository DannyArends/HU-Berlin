# Analysis of missing hetrozygosity in the Mega Muga array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))
plot.chr <- function(chrinfo){
  cnt <- 1; invisible(apply(chrinfo,1,function(x){ lines(c(cnt, cnt), c(0, x["Length"]), type="l", col="black", lty=1, lwd=2); cnt <<- cnt + 1 }))
}

setwd("E:/Mouse/DNA/DiversityArray/")
chrinfo    <- read.table("Annotation/mouseChrInfo.txt", colClasses=c("character", "numeric"), header=TRUE)                          # Chromosome information
max.length <- max(chrinfo[,"Length"])

setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                     # Read in the data from the Mega Muga
genotypes <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                                      # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),]                                               # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]
P  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]

setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                     # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes <- genotypes[,F2]

he <- apply(genotypes, 1, function(x){ sum(x == "H", na.rm=TRUE) })                                                                 # Calculate the number of heterozygous alleles
ho <- apply(genotypes, 1, function(x){ sum(x != "H", na.rm=TRUE) })                                                                 # Calculate the number of homozygous alleles

window.size <- 750000
step.size <- 750000/2
nperms <- 1000

bins <- NULL
for(chr in chromosomes){
  chr.length <- chrinfo[chrinfo[,"Chr"] == chr, "Length"]                                                                           # Length of the chromosome
  chr.bins <- cbind(Start = seq(1, chr.length, step.size), Stop = seq(1, chr.length, step.size) + window.size)                      # Create our bins
  smap <- map[map[,"Chr"] == chr, ]                                                                                                 # Create a copy of the map (only this chromosome)
  
  stats <- t(apply(chr.bins, 1, function(x){
    chr.markers <- rownames(smap)[as.numeric(smap[,"Mb_NCBI38"]) > x["Start"] & as.numeric(smap[,"Mb_NCBI38"]) < x["Stop"]]         # Which markers are in my bin
    return(c(length(chr.markers), sum(he[chr.markers]) / sum(ho[chr.markers])))                                                     # Calculate HE / HO ratio
  }))
  chr.bins <- cbind(chr.bins, nSNPs = stats[, 1], Score = stats[, 2])                                                               # Remember the scores

  for(perm in 1:nperms){                                                                                                            # permute this chromosome
    permstats <- t(apply(chr.bins, 1, function(x){
      if(as.numeric(x["nSNPs"]) != 0){                                                                                              # If there are SNPs in the bin
        chr.markers <- rownames(smap)[sample(length(rownames(smap)), as.numeric(x["nSNPs"]))]                                       # Draw N SNPs at random from this chromosome
        return(c(length(chr.markers), sum(he[chr.markers]) / sum(ho[chr.markers])))                                                 # Calculate HE / HO ratio
      }else{
        return(c(0, NaN))                                                                                                           # No SNPs, just return no score
      }
    }))
    chr.bins <- cbind(chr.bins, permstats[, 2])                                                                                     # Add the scores observed during permutation
  }
  cat(paste0("Done with chromosome ",chr,"/",length(chromosomes),"\n"))
  chr.bins <- cbind(Chr = chr, chr.bins)
  bins <- rbind(bins, chr.bins)
}
colnames(bins)[which(colnames(bins) == "")] <- paste0("P" , 1:length(which(colnames(bins) == "")))
permutationmatrix <- apply(bins[,grep("^P", colnames(bins))],2,as.numeric)

bin.thresholds <- t(apply(permutationmatrix,1,function(x){
  return(c(mean(x) + 5 * sd(x), mean(x) - 5 * sd(x)))
}))

ratio.above <- which(bins[,"Score"] > bin.thresholds[,1])
ratio.below <- which(bins[,"Score"] < bin.thresholds[,2])

regions.above <- apply(bins[ratio.above,1:3], 1, paste0, collapse=":")
regions.below <- apply(bins[ratio.below,1:3], 1, paste0, collapse=":")

bm.above <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "mgi_id","mgi_symbol", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region", "biotype"), values = list(regions.above,"protein_coding"), mart = bio.mart)
bm.below <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "mgi_id","mgi_symbol", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region", "biotype"), values = list(regions.below,"protein_coding"), mart = bio.mart)

cat("Found",length(unique(bm.above[,"ensembl_gene_id"])), "genes in hyper-heterozygous regions\n")
cat("Found",length(unique(bm.below[,"ensembl_gene_id"])), "genes in hypo-heterozygous regions\n")

cat(unique(bm.above[,"ensembl_gene_id"]),sep="\n")
cat(unique(bm.below[,"ensembl_gene_id"]),sep="\n")

### HWE

ps <- apply(genotypes, 1, function(x){
  obsAA <- sum(x == "A", na.rm=TRUE) ; obsH <- sum(x == "H", na.rm=TRUE); obsBB <- sum(x == "B", na.rm=TRUE)
  P <- (2 * obsAA + obsH) / (2 * ( obsAA + obsH + obsBB))
  Q <- (1 - P)
  n <- sum(!is.na(genotypes[1,]))
  expAA <- (P^2) * n; expH <- 2 * P * Q * n; expBB <- (Q^2) * n
  chisq <- ((obsAA - expAA)^2 / expAA) + ((obsH - expH)^2 / expH) + ((obsBB - expBB)^2 / expBB)
  pchisq(chisq, 1, lower=FALSE)
})





boxplot(t(permutationmatrix))
points(bins[,"Score"], col="green")

up2sd <- apply(permutationmatrix,1,mean) + (2 * apply(permutationmatrix, 1, sd)) ; up3sd <- apply(permutationmatrix,1,mean) + (5 * apply(permutationmatrix, 1, sd))
dw2sd <- apply(permutationmatrix,1,mean) - (2 * apply(permutationmatrix, 1, sd)) ; dw3sd <- apply(permutationmatrix,1,mean) - (5 * apply(permutationmatrix, 1, sd))
dw2sd[dw2sd < 0] <- 0 ; dw3sd[dw3sd < 0] <- 0 ;

co <- rep("black", nrow(bins))
#co[which(bins[,"Score"] > up2sd)] <- "blue" ; co[which(bins[,"Score"] < dw2sd)] <- "orange"
co[which(bins[,"Score"] > up3sd)] <- "green" ; co[which(bins[,"Score"] < dw3sd)] <- "red"

plot(c(1,nrow(bins)), c(0,3), t = 'n')
points(up3sd, t = 'l', col="orange",lwd=0.5,lty=1)
points(dw3sd, t = 'l', col="orange",lwd=0.5,lty=1)
points(bins[,"Score"], col=co, pch=20,cex=0.8)

plot(y=c(0, max.length), x=c(1,nrow(chrinfo)), t='n', yaxt="n", xlab="", ylab="", xaxt="n")
abline(h=seq(0, max.length, 10000000), col = "lightgray", lty = "dotted")

plot.chr(chrinfo)

cnt <- 1
aa <- apply(bins,1,function(x){
  xloc <- match(as.character(x["Chr"]), chromosomes); yloc <- as.numeric(x["Start"])
  if(co[cnt] != "black") points(x=xloc, y=yloc, pch="-", col=co[cnt], cex=2.5)
  cnt <<- cnt + 1
})

cnt <- 1
aa <- lapply(-log10(ps), function(x){
  xloc <- match(as.character(map[cnt,"Chr"]), chromosomes); yloc <- as.numeric(map[cnt,"Mb_NCBI38"])
  if(!is.nan(x) && x > 5) points(x=xloc+0.1, y=yloc, pch="-", cex=1.0)
  cnt <<- cnt + 1
})

axis(1,chrinfo[,1], at=c(1:nrow(chrinfo)), las = 1, cex.axis = 1.5)
axis(2, seq(0, max.length, 10000000)/1000000, at = seq(0,  max.length, 10000000), cex.axis = 1.2, las = 1)

#### OLD

regions.red <- apply(bins[which(co=="red"),1:3], 1, paste0, collapse=":")
regions.green <- apply(bins[which(co=="green"),1:3], 1, paste0, collapse=":")

res.red <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "mgi_id","mgi_symbol", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region", "biotype"), values = list(regions.red,"protein_coding"), mart = bio.mart)
res.green <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "mgi_id","mgi_symbol", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region", "biotype"), values = list(regions.green,"protein_coding"), mart = bio.mart)

interactions <- read.csv("Additional/BIOGRID-ALL-3.2.121.tab2.txt", sep="\t", header=TRUE)
interactions <- interactions[which(interactions[,"Organism.Interactor.A"] == 10090),c("Entrez.Gene.Interactor.A", "Entrez.Gene.Interactor.B")]

ints <- which(interactions[,"Entrez.Gene.Interactor.A"] %in% res.green[,"entrezgene"] & interactions[,"Entrez.Gene.Interactor.B"] %in% res.green[,"entrezgene"] & interactions[,"Entrez.Gene.Interactor.A"] != interactions[,"Entrez.Gene.Interactor.B"])
interactions <- interactions[ints,]

entreztoloc <- function(entrezgene = 17060){ res.green[which(res.green[,"entrezgene"] == entrezgene),c("chromosome_name", "start_position", "end_position")] }

circlelocations <- function(nt){
  medpoints <- matrix(nrow = nt, ncol = 2) ; phi <- seq(0, 2 * pi, length = (nt + 1)) ; complex.circle <- complex(modulus = 1, argument = phi)
  for(j in 1:nt){ medpoints[j, ] <- c(Im(complex.circle[j]), Re(complex.circle[j])) }
  return(medpoints)
}

drawspline <- function(cn1, cn2, lwd = 1,col="blue",...){
  x <- cbind(cn1[1],0,cn2[1]) ; y <- cbind(cn1[2],0,cn2[2]); r <- xspline(x, y, lty=1, shape=1, lwd=lwd, border=col,...)
}

gapsize <- 30
chromosome <- 5

clocs <- circlelocations(ceiling(sum(chrinfo[,"Length"]) / 1000000) + nrow(chrinfo) * gapsize)

toCircle <- function(chrinfo, chr, mb = 5, gapsize = 10){
  ix <- which(chrinfo[,"Chr"] == chr)-1; t(clocs[(sum(chrinfo[0:ix,"Length"])/1000000) + mb + gapsize * ix,])
}

plot(c(-1.1, 1.1), c(-1.1, 1.1), type = "n", axes = FALSE, xlab = "", ylab = "")
l <- 1
for(chr in chromosomes[1:21]){
  nl <- l + round(chrinfo[chrinfo[,"Chr"] == chr,"Length"] / 1000000)
  lines(clocs[l:nl,],cex=0.02)
  l <- nl + gapsize
}

cnt <- 1
aa <- apply(bins,1,function(x){
  loc <- toCircle(chrinfo, x["Chr"], as.numeric(x["Start"])/1000000, gapsize)
  if(co[cnt] != "black") lines(rbind(0.98*loc,1.02*loc), col=co[cnt], cex=1)
  cnt <<- cnt + 1
})

for(i in 1:nrow(interactions)){
  loc1 <- entreztoloc(interactions[i,"Entrez.Gene.Interactor.A"])
  loc2 <- entreztoloc(interactions[i,"Entrez.Gene.Interactor.B"])
  #if(loc1["chromosome_name"] == chromosome | loc2["chromosome_name"] == chromosome){
    cloc1 <- toCircle(chrinfo, as.character(loc1["chromosome_name"]), as.numeric(loc1["start_position"])/1000000, gapsize)
    cloc2 <- toCircle(chrinfo, as.character(loc2["chromosome_name"]), as.numeric(loc2["start_position"])/1000000, gapsize)
    drawspline(cloc1, cloc2, lwd = 1, col="gray")
  #}
}


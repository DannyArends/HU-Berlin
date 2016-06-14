# Analysis of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2016
# first written Jun, 2016

### Load in SNP data, sample and SNP information files

setwd("E:/Goat/DNA/SihamAnalysis")

numsnpdata <- read.csv("filtered_snps_numeric_NO_DN2.txt", sep="\t", check.names=FALSE)
snpinfo <- read.csv("merged_snp_info.txt", sep="\t", check.names=FALSE)
samples <- read.csv("merged_samples_NO_DN2.txt", sep="\t", check.names=FALSE)

# Remove non-seggregating markers ( 1 is left)

toR <- names(which(unlist(lapply(apply(numsnpdata,1,table),length)) == 1))
snpinfo <- snpinfo[-which(rownames(snpinfo) == toR),]
numsnpdata <- numsnpdata[-which(rownames(numsnpdata) == toR),]

pheNames <- c("Averagemilk", "Weight", "Withersheight", "Rumpheight", "Bodylength", "Sternumheight", 
              "Bodydepth", "Bicoastaldiameter", "Earlength", "RumpWidth", "HeadWidth", "Rumplength", 
              "Headlength", "Heartgirth", "Cannonbone", "Muzzlediameter")

phenotypes <- samples[colnames(numsnpdata),pheNames]
breed  <- as.factor(samples[colnames(numsnpdata),"Breed"])
location  <- as.factor(samples[colnames(numsnpdata),"locationShort"])

setwd("E:/Goat/DNA/SihamQTL")

### Create supplement table 1

breedeffects <- NULL
meansandsds <- NULL
for(x in 1:ncol(phenotypes)){
  breedeffects <- rbind(breedeffects, c(colnames(phenotypes)[x], anova(lm(phenotypes[,x] ~ breed))[[5]][1]))
  tagg <- round(c(mean(phenotypes[breed == "Tagg",x],na.rm=TRUE), sd(phenotypes[breed == "Tagg",x],na.rm=TRUE)),3)
  nu <- round(c(mean(phenotypes[breed == "Nu",x],na.rm=TRUE), sd(phenotypes[breed == "Nu",x],na.rm=TRUE)),3)
  ni <- round(c(mean(phenotypes[breed == "Ni",x],na.rm=TRUE), sd(phenotypes[breed == "Ni",x],na.rm=TRUE)),3)
  dese <- round(c(mean(phenotypes[breed == "Dese",x],na.rm=TRUE), sd(phenotypes[breed == "Dese",x],na.rm=TRUE)),3)
  meansandsds <- rbind(meansandsds, c(tagg, ni, nu, dese))
}
breedeffects <- cbind(breedeffects, p.adjust(breedeffects[,2], "BH"), meansandsds)
colnames(breedeffects) <- c("Phenotype", "Pvalue", "Padjusted", "mean(Tagg)", "sd(Tagg)", "mean(Ni)", "sd(Ni)", "mean(Nu)", "sd(Nu)", "mean(Dese)", "sd(Dese)")
write.table(breedeffects, "EffectOfBreed.txt", sep = "\t",row.names=FALSE,quote=FALSE)

### Fancy phenotype plot (Figure 1)

rot <- seq(1, 360, 360/length(pheNames))                                                        # Rotational axis
names(rot) <- pheNames                                                                            # Phenotype names displayed for each axis
cp <- seq(0, 2 * pi, length = 360)                                                                # Circle points in degrees
coords <- t(rbind( sin(cp), cos(cp)))                                                             # Circle coordinates
colz <- c(rgb(1,0,0,0.3), rgb(0,1,0,0.3), rgb(0,0,1,0.3), rgb(1,0,1,0.3)) 
colzNT <- c(rgb(1,0,0), rgb(0,1,0), rgb(0,0,1), rgb(1,0,1)) 
breeds <-  as.character(unique(breed))
cnt <- 1

op <- par(mfrow=c(2, 2), mai = c(0.05, 0.1, 0.2, 0.1))                                            # 6 sublines ( so we use a 3x2 plot )

for(b in breeds){
  plot(c(-1.3, 1.3),c(-1.3, 1.3), t = 'n', ylab="", xlab="",xaxt='n',yaxt='n', main=b)  # Create the plot window
  points(coords, pch = 18, cex = 0.4)                                                             # Add the circle
  for(x in 1:length(rot)){
    points(c(0, coords[rot[x],1]), c(0,coords[rot[x],2]), t='l', lwd=0.5, lty=3)                  # The different phenotype axis
    text(1.2 * coords[rot[x],1], 1.1 * coords[rot[x],2], names(rot)[x], cex=0.8)                  # The name on the outside of the circle
  }
  scoords <- NULL
  mincoords <- NULL
  maxcoords <- NULL
  for(phe in names(rot)){
    minP <- min(as.numeric(phenotypes[,phe]),na.rm=TRUE)
    maxP <- max(as.numeric(phenotypes[,phe]),na.rm=TRUE)
    meanP <- mean(as.numeric(phenotypes[breed == b, phe]),na.rm=TRUE)
    sdP <- sd(as.numeric(phenotypes[breed == b, phe]),na.rm=TRUE)
    phescale <- ((meanP - minP) / (maxP - minP))                                                # Scale
    pheMin <- (((meanP-sdP) - minP) / (maxP - minP))                                                # Scale
    pheMax <- (((meanP+sdP) - minP) / (maxP - minP))                                                # Scale
    coord <- c(phescale * coords[rot[phe],1], phescale * coords[rot[phe],2])                    # X and Y coordinate on the scaled circle
    coordMin <- c(pheMin * coords[rot[phe],1], pheMin * coords[rot[phe],2])                    # X and Y coordinate on the scaled circle
    coordMax <- c(pheMax * coords[rot[phe],1], pheMax * coords[rot[phe],2])                    # X and Y coordinate on the scaled circle
 #   points(rbind(t(coordMin),t(coordMax)), pch=19,cex=0.6, col = colzNT[cnt], t ='l', lwd=2)                                         # plot the dot
    points(t(coord), pch=19,cex=0.8, col = colzNT[cnt])                                         # plot the dot
    
    mincoords <- rbind(mincoords, coordMin)    
    maxcoords <- rbind(maxcoords, coordMax)    
    
    scoords <- rbind(scoords, coord)                                                            # remember the coordinate

  }
  rownames(scoords) <- names(rot)                                                               # Set rownames on the polygon edges
  polygon(scoords, col = colz[cnt], pch=18,fillOddEven = TRUE)                                  # Connect all the dots using a polygon
  polygon(maxcoords, col = colz[cnt], pch=18,fillOddEven = TRUE)                                  # Connect all the dots using a polygon
  polygon(mincoords, col = "white", pch=18,fillOddEven = TRUE)                                  # Connect all the dots using a polygon
  cnt <- cnt + 1
}

### Correlation plot

op <- par(mai = c(2,2,1,1))
corP <- cor(phenotypes, use="pair")
image(abs(corP), xaxt='n', yaxt='n', breaks = seq(0,1,0.1), col=gray(seq(1.0, 0.1, -0.1)))
axis(1, at = seq(0,1,1/(ncol(corP)-1)), rownames(corP),las=2,cex.axis=0.8)
axis(2, at = seq(0,1,1/(ncol(corP)-1)), colnames(corP),las=2,cex.axis=0.8)
abline(h = seq(0,1,1/(ncol(corP)-1)) + 0.5 * 1/(ncol(corP)-1), col="white")
abline(v = seq(0,1,1/(ncol(corP)-1)) + 0.5 * 1/(ncol(corP)-1), col="white")
for(x in 1:ncol(corP)){
  for(y in 1:nrow(corP)){
    col <- "white"; if(abs(corP[x,y]) < 0.4) col <- "black"
    text(1/(ncol(corP)-1) * (x-1), 1/(ncol(corP)-1) * (y-1), round(corP[x,y],2), col=col)
  }
}
box()

### QTL analysis


## Breed in the model
for(p in 1:ncol(phenotypes)){
  pvalues <- NULL
  for(x in 1:nrow(numsnpdata)){
    pvals <- anova(lm(phenotypes[,p] ~ breed + as.numeric(numsnpdata[x,])))[[5]]
    if(length(pvals) == 3) {
      pvalues <- rbind(pvalues, c(rownames(numsnpdata)[x], pvals))
    }else{
      cat("Unable to fit model at marker ",x," :", rownames(numsnpdata)[x], "\n")
    }
  }
  op <- par(mfrow=c(1,2))
  plot(-log10(as.numeric(pvalues[,3])))
  plot(-log10(as.numeric(pvalues[,2])))
  cat(colnames(phenotypes)[p], " = ", length(which(p.adjust(as.numeric(pvalues[,3]), "BH") < 0.05)), " marker\n")
  write.table(pvalues, file=paste0(colnames(phenotypes)[p],"_qtl.txt"), sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE)
}

for(p in 1:ncol(phenotypes)){
  results <- read.table(paste0(colnames(phenotypes)[p],"_qtl.txt"), sep="\t")
  op <- par(mfrow=c(1,2))
  plot(-log10(as.numeric(results[,3])))
  plot(-log10(as.numeric(results[,2])))
  cat(colnames(phenotypes)[p], " = ", max(-log10(as.numeric(results[,3]))), " marker\n")
}

## Breed corrected first, then basic QTL mapping

phenotypesN <- phenotypes
for(p in 1:ncol(phenotypes)){
  res <- lm(phenotypes[,p] ~ breed)
  phenotypesN[,p] <- NA
  phenotypesN[as.numeric(names(res$residuals)), p] <- res$residuals
}

## Breed in the model
for(p in 1:ncol(phenotypesN)){
  pvalues <- NULL
  for(x in 1:nrow(numsnpdata)){
    pvals <- anova(lm(phenotypesN[,p] ~ breed + as.numeric(numsnpdata[x,])))[[5]]
    if(length(pvals) == 3) {
      pvalues <- rbind(pvalues, c(rownames(numsnpdata)[x], pvals))
    }else{
      cat("Unable to fit model at marker ",x," :", rownames(numsnpdata)[x], "\n")
    }
  }
  op <- par(mfrow=c(1,2))
  plot(-log10(as.numeric(pvalues[,3])))
  plot(-log10(as.numeric(pvalues[,2])))
  cat(colnames(phenotypesN)[p], " = ", length(which(p.adjust(as.numeric(pvalues[,3]), "BH") < 0.05)), " marker\n")
  write.table(pvalues, file=paste0(colnames(phenotypesN)[p],"_adj_qtl.txt"), sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE)
}
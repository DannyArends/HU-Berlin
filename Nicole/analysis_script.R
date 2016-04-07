#
# Analysis of Nicole's 'complicated QTL data' using multiple QTL mapping
# copyright (c) 2016-2020 - Danny Arends and Nicole Hallahan
# last modified Apr, 2016
# first written Apr, 2016
#

##### Step 1: Data loading an processing into the correct format

setwd("D:/Collegues/Nicole")                                                # Change this to the folder where data.txt is stored
mdata <- read.table("data.txt", sep="\t")                                   # Loads all the data into R

map <- mdata[2, c(36:87)] ; map <- map[-c(1:4)]                             # split out the map, (throw away the first 4 markers)
locs <- as.numeric(as.character(unlist(map)))                               # locations on the map

genotypes <- mdata[3:203, c(36:87)] ; genotypes <- genotypes[,-c(1:4)]      # split out the genotypes, (throw away the first 4 markers)
genotypes[169,"V84"] <- "ND"                                                # fix ?! a single genotype value

numGeno <- apply(genotypes,2,function(x){as.numeric(as.factor(x))})         # genotypes as numbers (not factors)

# Subset the phenotypes, so that we have the most interesting ones.
phenotypes <- apply(mdata[3:203, c(11, 22, 25, 26, 29:34)], 2, function(x){as.numeric(x)})
colnames(phenotypes) <- c("BW week 13", "BG week 13", "BG AUC wk 11", "Ins AUC wk 11", 
                          "Final plasma TG", "Glycerol", "FFA", "Final plasma CHOL", 
                          "Liver mg TG per mg prot", "Liver mg CHOL per mg prot")

phenotypes <- phenotypes[,1:2]

##### Step 2: Definition of different functions, we use to do QTL mapping

# Function that extracts the last LOD score from a multiple QTL model
aslodscores <- function(qtldata) {
  -log10(unlist(lapply(qtldata[["models"]], function(model){
    # cat("model-length:", (length(x[[5]])-1), "\n"); 
    x <- anova(model)
    return(x[[5]][(length(x[[5]])-1)]);
  })))
}

# Function that extracts the variance explained from a multiple QTL model
asvarexplained <- function(qtldata) {
  lapply(qtldata[["models"]], function(model){
    # cat("model-length:", (length(x[[5]])-1), "\n"); 
    x <- anova(model)
    100* (x[[2]] / sum(x[[2]]))
  })
}

# Function to map a QTL using the model: Y = Cmark + marker + err
# Set cmark1 and scan for a second QTL using a multiple QTL model, drop the marker covariate when it is nearby
QTLmappingOne <- function(genotypes, phenotypes, locs, cpheno =  2, cmark = 10) {
  nearby <- which(locs > (locs[cmark] - 10) & locs < (locs[cmark] + 10))
  mcnt <- 0
  models <- apply(genotypes, 2, function(x){
    marker <- as.factor(as.character(x))
    m1 <- as.factor(as.character(genotypes[,cmark]))
    mcnt <<- mcnt + 1
    if(mcnt %in% nearby){
      return(lm(phenotypes[,cpheno] ~ marker))
    }else{
      return(lm(phenotypes[,cpheno] ~ m1 + marker))
    }
  })
  colz <- rep("orange", nrow(genotypes))
  colz[nearby] <- rep("blue", length(nearby))
  colz[cmark] <- "green"
  invisible(list(models = models, locs = locs, nearby = nearby, cmark = cmark, colorcodes = colz))
}

# Draw the LOD scores of the single QTL model, locs = map location of the markers, lodscores contains the LOD scores
# Warning: This function does not setup a plot window, it assumes a window is already there
plotBaseLine <- function(locs, lodscores){
  points(locs, lodscores, t = 'l')
  points(locs, lodscores, t = 'p', pch=19, cex=0.5)
  points(locs, rep(0, length(lodscores)), pch="|")
}

# Draw the output of a multiple QTL model (model: Y = Cmark + marker + err)
# Warning: This function does not setup a plot window, it assumes a window is already there
plotmodel1 <- function(qtldata) {
  lods <- aslodscores(qtldata)
  abline(v = qtldata[["locs"]][ min(qtldata[["nearby"]]) ]-1, col="gray", lty=2, lwd=0.4)
  abline(v = qtldata[["locs"]][ max(qtldata[["nearby"]]) ]+1, col="gray", lty=2, lwd=0.4)
  points(qtldata[["locs"]], lods, t = 'l', col='orange', lwd = 2)
  points(qtldata[["locs"]], lods, t = 'p', pch=19, cex=0.8, col=qtldata[["colorcodes"]])
  legend("topright", c("Single QTL model", "Multiple QTL model", "nearby", "cofactor"), col=c("black", "orange", "blue", "green"), lwd=c(2, 2, NA, NA), pch=c(19,19,19,19), cex=0.8)
}

##### Step 3: do the actual single QTL mapping and the  multiple QTL model (model: Y = Cmark + marker + err)

#for(phe in 1:ncol(phenotypes)){      # Uncomment this to make all the plots for all the phenotypes
phe <- 2                              # Set it here if you want to analyze a single phenotype to make a specific plot

  models <- apply(genotypes, 2, function(x){
    marker <- as.factor(as.character(x))
    return(anova(lm(phenotypes[,phe] ~ marker)))
  })

  lods <- -log10(unlist(lapply(models, function(x){ x[[5]][1]; })))
  vars <- lapply(models, function(x){  100* (x[[2]] / sum(x[[2]])) })
  for(x in 1:length(vars)){
    cat(colnames(phenotypes)[phe], locs[x], " variance explained (single QTL)", vars[[x]], "\n")
  }
  # Print some information for the table
  #cat(phe, "", max(lods), "", locs[which.max(lods)], locs[ which(lods > (max(lods)-2)) ], "\n")

  op <- par(mfrow=c(1,1))
  plot(c(80, 160), c(0, 10), t = 'n', xlab="Location (mBp)", ylab = "LOD score -log10(pvalue)")
  plotBaseLine(locs, lods)
  
  # --> If you comment the next line there will be a plot in R, otherwise on the HDD, also comment the dev.off() call below
  png(paste0(colnames(phenotypes)[phe], ".png"), width=1280, height=800)
  op <- par(mfrow=c(2, 2))
  for(cmark in c(5, 25, 32, 48)) { # Update here to move the markers
    plot(c(80, 160), c(0, 10), t = 'n', xlab="Location (mBp)", ylab = "LOD score -log10(pvalue)")
    plotBaseLine(locs, lods)
    qtldata <- QTLmappingOne(genotypes, phenotypes, locs, cpheno = phe, cmark = cmark)
    plotmodel1(qtldata)
    vExplained <- asvarexplained(qtldata)
    for(y in 1:length(vExplained)){
      cat(colnames(phenotypes)[phe], locs[cmark], locs[y], "variance explained", vExplained[[y]], "\n")
    }
  }
  # --> the next line shouldshould match the png line.
  dev.off()                                                                   
  cat("Done", colnames(phenotypes)[phe], "\n")
#}

# Correlation plot of genotypes against eachother
corM <- cor(numGeno, use="pair")

image(1:nrow(corM), 1:nrow(corM), corM, xaxt='n', yaxt='n', xlab="", ylab="")     # Plot the correlation matrix
axis(1, at=1:nrow(corM), locs, las=2, cex=0.8)                                    # Add locations to the x-axis
axis(2, at=1:nrow(corM), locs, las=2, cex=0.8)                                    # Add locations to the y-axis
box()

##### Step 4: function for the multiple QTL model (model: Y = Cmark1 + Cmark2 + marker + err)

# Code, using 2 markers as covariates, and mapping the genetic map using a 3 QTL model
QTLmappingTwo <- function(genotypes, phenotypes, locs, cpheno =  2, cmark1 = 5, cmark2 = 40) {
  nearby1 <- which(locs > (locs[cmark1] - 10) & locs < (locs[cmark1] + 7.3))    # UPDATE: To play with the area around marker1
  nearby2 <- which(locs > (locs[cmark2] - 13) & locs < (locs[cmark2] + 20))     # UPDATE: To play with the area around marker2
  nearby3 <- which(locs > (locs[48] - 5) & locs < (locs[48] + 10))     # UPDATE: To play with the area around marker2
  mcnt <- 0
  models <- apply(genotypes, 2, function(x){
    marker <- as.factor(as.character(x))
    m1 <- as.factor(as.character(genotypes[,cmark1]))
    m2 <- as.factor(as.character(genotypes[,cmark2]))
    mcnt <<- mcnt + 1
    if(mcnt %in% nearby1) {
      return(lm(phenotypes[,cpheno] ~ m2 + marker))
    } else if(mcnt %in% nearby2) {
      return(lm(phenotypes[,cpheno] ~ m1 + marker))
    } else {
      return(lm(phenotypes[,cpheno] ~ m2 + m1 + marker))
    }
  })
  colz <- rep("orange", length(lods))
  colz[nearby1] <- rep("blue", length(nearby1))
  colz[nearby2] <- rep("purple", length(nearby2))
  colz[cmark1] <- "green"
  colz[cmark2] <- "red"
  invisible(list(models = models, locs = locs, nearby1 = nearby1, nearby2 = nearby2, cmark1 = cmark1, cmark2 = cmark2, colorcodes = colz))
}

# Draw the output of the multiple QTL model (model: Y = Cmark1 + Cmark2 + marker + err)
plotmodel2 <- function(qtldata) {
  lods <- aslodscores(qtldata)
  abline(v = qtldata[["locs"]][ min(qtldata[["nearby1"]]) ]-1, col="gray", lty=2, lwd=0.4)
  abline(v = qtldata[["locs"]][ max(qtldata[["nearby1"]]) ]+1, col="gray", lty=2, lwd=0.4)
  abline(v = qtldata[["locs"]][ min(qtldata[["nearby2"]]) ]-1, col="gray", lty=2, lwd=0.4)
  abline(v = qtldata[["locs"]][ max(qtldata[["nearby2"]]) ]+1, col="gray", lty=2, lwd=0.4)
  points(qtldata[["locs"]], lods, t = 'l', col='orange', lwd = 2)
  points(qtldata[["locs"]], lods, t = 'p', pch=19, cex=0.8, col=qtldata[["colorcodes"]])
  legend("topright", c("Single QTL model", "Multiple QTL model", "nearby", "cofactor"), col=c("black", "orange", "blue", "green"), lwd=c(2, 2, NA, NA), pch=c(19,19,19,19), cex=0.8)
  abline(h = -log10(0.05/50), lty=2, col="green")
}

phe <- 2 # Set the phenotype to 2, and remap the single marker model for the phenotype
models <- apply(genotypes, 2, function(x){
  marker <- as.factor(as.character(x))
  return(anova(lm(phenotypes[,phe] ~ marker)))
})

# Lod scores from the single marker model and QTL mapping using a 3 QTL model
lods <- -log10(unlist(lapply(models, function(x){ x[[5]][1]; })))
qtldata <- QTLmappingTwo(genotypes, phenotypes, locs, cpheno = phe, cmark1 = 4, cmark2 = 37)  # Update to put your covariates at a differetn marker

# Create the plot
#png("multiple QTL model_BG.png", width=1280, height=800)  # If you need to write it out to disk
  plot(c(80, 160), c(0, 10), t = 'n', xlab="Location (mBp)", ylab = "LOD score -log10(pvalue)")
  plotBaseLine(locs, lods)
  plotmodel2(qtldata)
#dev.off()









###---- OLD Code down here, beware !
plot(c(80,154), c(0,10))
points(as.numeric(as.character(unlist(map))), -log10(unlist(lapply(models, function(x){
  x[[5]][2]
}))), t = 'l')

points(as.numeric(as.character(unlist(map))), -log10(unlist(lapply(models, function(x){
  x[[5]][1]
}))), t = 'l', col='red')

models[[25]] # 25 is in the middle of the QTL region at 110 cm
models[[4]] # 25 is in the middle of the QTL region at 110 cm
100* (models[[25]][[2]] / sum(models[[25]][[2]]))

# Marker 4 is near the frontal QTL peek at 83.6 cM
models <- apply(genotypes, 2, function(x){
  marker <- as.factor(as.character(x))
  m1 <- as.factor(as.character(genotypes[,4]))
  if(length(levels(marker)) > 1){
    return(anova(lm(phenotypes[,phe] ~ m1 + marker)))
  }
})


plot(c(80,154), c(0,10))

points(as.numeric(as.character(unlist(map))), -log10(unlist(lapply(models, function(x){
  x[[5]][2]
}))), t = 'l')

#Region distal: from 95.40 to 120

points(as.numeric(as.character(unlist(map))), -log10(unlist(lapply(models, function(x){
  x[[5]][1]
}))), t = 'l', col='red')

# 25 is in the middle of the QTL region at 110 cm
100* (models[[25]][[2]] / sum(models[[25]][[2]]))


# Marker 4 and 35
models <- apply(genotypes, 2, function(x){
  marker <- as.factor(as.character(x))
  m1 <- as.factor(as.character(genotypes[,4]))
  m2 <- as.factor(as.character(genotypes[,35]))
  if(length(levels(marker)) > 1){
    return(anova(lm(phenotypes[,phe] ~ m1 + m2 + m1:m2 + marker)))
  }
})

pvals <- NULL
for(x in 1:length(map)){
  if(x == 4) {
    pvals <- c(pvals, models[[x]][[5]][1])
  }else if(x == 35){
    pvals <- c(pvals, models[[x]][[5]][2])
  }else{
    pvals <- c(pvals, models[[x]][[5]][3])
  }  
}
plot(-log10(pvals), t = 'l')

AIC(models[[25]])

plot(c(80,154), c(0,10))

points(as.numeric(as.character(unlist(map))), -log10(unlist(lapply(models, function(x){
  x[[5]][3]
}))), t = 'l')

points(as.numeric(as.character(unlist(map))), -log10(unlist(lapply(models, function(x){
  x[[5]][2]
}))), t = 'l', col='blue')

points(as.numeric(as.character(unlist(map))), -log10(unlist(lapply(models, function(x){
  x[[5]][1]
}))), t = 'l', col='red')

# 25 is in the middle of the QTL region at 110 cm
100* (models[[25]][[2]] / sum(models[[25]][[2]]))


m1 <- as.factor(as.character(genotypes[,8]))   #ND = high BG, NN low BG, 94.9
m2 <- as.factor(as.character(genotypes[,35]))  #ND = high BG, NN low BG, 121
ordering <- which(m1 != m2)

genotypes <- genotypes[ordering, ]
m3 <- as.factor(as.character(genotypes[,20]))   #ND = high BG, NN low BG, 105.27
m4 <- as.factor(as.character(genotypes[,37]))   #ND = high BG, NN low BG, 127.7

ordering <- which(m3 != m4)

numgeno <- apply(genotypes[ordering,], 2,function(x){ return(as.numeric(as.factor(as.character(x)))) })
order2 <- sort(phenotypes[ordering,2], index.return=TRUE)$ix
image(1:ncol(numgeno),1:nrow(numgeno), t(numgeno[order2,]))
text(rep(2,length(phenotypes[ordering,2][order2])), 1:length(phenotypes[ordering,2][order2]), phenotypes[ordering,2][order2], cex=0.6)
text(1:ncol(genotypes), rep(2, ncol(genotypes)),  as.numeric(as.character(unlist(map))), cex=0.6,las=2)



genotypes[-which(m2!="NN"),]

high <- which(phenotypes[,2] > 16)
low <- which(phenotypes[,2] < 11)
medium <- which(!(1:length(phenotypes) %in% c(high,low)))


#Something Something check bodyweight

phenotypes <- apply(mdata[3:203, c(11,22)], 2, function(x){as.numeric(x)})
phenotypes <- apply(mdata[3:203, c(1,2,3,4,5,6,7,8,9,10)], 2, function(x){as.numeric(x)})
#phenotypes <- apply(mdata[3:203, c(12,13,14,15,16,17,18,19,20,21)], 2, function(x){as.numeric(x)})

for( phe in 2:10){
  models <- apply(genotypes, 2, function(x){
    marker <- as.factor(as.character(x))
    if(length(levels(marker)) > 1){
      return(anova(lm((phenotypes[,phe] - phenotypes[,(phe-1)]) ~ marker)))
    }
  })

  lods <- -log10(unlist(lapply(models, function(x){
    x[[5]][1]
  })))
  op <- par(mfrow=c(1,2))
  plot(as.numeric(as.character(unlist(map))), lods, t = 'l')
  plot((phenotypes[,phe] - phenotypes[,(phe-1)]) ~ genotypes[,which.max(lods)])
  cat(phe, max(lods), "\n")
  readline()
}





# No covariate
models <- apply(genotypes, 2, function(x){
  marker <- as.factor(as.character(x))
  if(length(levels(marker)) > 1){
    return(anova(lm(phenotypes[,phe] ~ marker)))
  }
})

pvals <- NULL
for(x in 1:length(map)){
#  if(x == 4) {
#    pvals <- c(pvals, models[[x]][[5]][1])
#  }else if(x == 35){
#    pvals <- c(pvals, models[[x]][[5]][2])
#  }else{
    pvals <- c(pvals, models[[x]][[5]][1])
#  }  
}
op <- par(mfrow=c(4,1))
plot(-log10(pvals), t = 'l')


# Marker 4
models <- apply(genotypes, 2, function(x){
  m1 <- as.factor(as.character(genotypes[,4]))
  marker <- as.factor(as.character(x))
  if(length(levels(marker)) > 1){
    return(anova(lm(phenotypes[,phe] ~ m1 + marker)))
  }
})

pvals <- NULL
for(x in 1:length(map)){
  if(x == 4) {
    pvals <- c(pvals, models[[x]][[5]][1])
#  }else if(x == 35){
#    pvals <- c(pvals, models[[x]][[5]][2])
  }else{
    pvals <- c(pvals, models[[x]][[5]][2])
  }  
}
plot(-log10(pvals), t = 'l')


# Marker 35
models <- apply(genotypes, 2, function(x){
  m1 <- as.factor(as.character(genotypes[,25]))
  marker <- as.factor(as.character(x))
  if(length(levels(marker)) > 1){
    return(anova(lm(phenotypes[,phe] ~ m1 + marker)))
  }
})

pvals <- NULL
for(x in 1:length(map)){
  #f(x == 4) {
   # pvals <- c(pvals, models[[x]][[5]][1])
#  }else 
  if(x == 25){
    pvals <- c(pvals, models[[x]][[5]][1])
  }else{
    pvals <- c(pvals, models[[x]][[5]][2])
  }  
}
plot(-log10(pvals), t = 'l')




# Marker 4 and 35
models <- apply(genotypes, 2, function(x){
  m1 <- as.factor(as.character(genotypes[,4]))
  m2 <- as.factor(as.character(genotypes[,35]))
  marker <- as.factor(as.character(x))
  if(length(levels(marker)) > 1){
    return(anova(lm(phenotypes[,phe] ~ m2 + m1 + marker)))
  }
})

pvals <- NULL
for(x in 1:length(map)){
  if(x == 4) {
    pvals <- c(pvals, models[[x]][[5]][2])
  }else if(x == 35){
    pvals <- c(pvals, models[[x]][[5]][1])
  }else{
    pvals <- c(pvals, models[[x]][[5]][3])
  }  
}
plot(-log10(pvals), t = 'l')


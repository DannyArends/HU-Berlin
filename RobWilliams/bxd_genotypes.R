#
# Analysis of BxD genotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

setwd("E:/Mouse/BxD")
alldata <- read.csv("genotypes.txt", sep = "\t", skip = 1, row.names=2, header=TRUE, na.strings=c("NA", "-3", ""))

physicalmap <- alldata[,c(3, 5, 7)]                                                     # extract the physical map
genotypes <- alldata[,-c(1:11)][,1:203]                                                 # extract the genotypes

genotypes["rs33672596","BXD146"] <- NA
genotypes["rs31572894","BXD93"]  <- NA
genotypes["rs31153518","BXD93"]  <- NA

toNumGeno <- function(genotypes){
  numgeno <- apply(genotypes, 2, function(x){ return(as.numeric(as.character(x))) })      # transform genotypes to numeric values
  rownames(numgeno) <- rownames(genotypes)                                                # Set the markernames
  return(numgeno)
}

cat("Starting with", nrow(genotypes), "of which", length(which(physicalmap[,3] == "Y")), "are selected by rob as markers\n")

#ii <- which(is.na(numgeno))[which(!(which(is.na(numgeno)) %in% which(is.na(genotypes))))]
#cat(ii, "\n")

calcRecombinations <- function(numgeno){
  recombinations <- vector("list", ncol(numgeno))
  for(i in 1:ncol(numgeno)){
    cgeno <- numgeno[1, i]
    recombination <- NULL
    for(m in 2:nrow(numgeno)){
      if(!is.na(numgeno[m, i]) && cgeno != numgeno[m, i]){
        if(physicalmap[rownames(numgeno)[m], "Chr"] == physicalmap[rownames(numgeno)[m-1], "Chr"]){
          recombination <- c(recombination, (m + (m-1)) / 2)
        }
        cgeno <- numgeno[m, i]
      }
    }
    recombinations[[i]] <- recombination
  }
  names(recombinations) <- colnames(genotypes)
  recombinationload <- unlist(lapply(recombinations, length))
  names(recombinationload) <- colnames(genotypes)
  return(list(recombinations = recombinations, recombinationload = recombinationload))
}

clean <- function(genotypes, recombinations){
  newgenotypes <- genotypes
  ll <- 0
  for(strain in colnames(genotypes)){
    ii <- which(diff(recombinations[[strain]]) < 5)
    for(x in ii){
      startloc <- (recombinations[[strain]][x] - 0.5)
      endloc <- (recombinations[[strain]][x+1] + 0.5)
      posdiff <-  max(physicalmap[startloc:endloc,2], na.rm=TRUE) - min(physicalmap[startloc:endloc,2], na.rm=TRUE)
      if(posdiff < 2){
        affected <- genotypes[startloc:endloc, strain]
        surrounding <- ceiling(length(affected) / 2)
        affectedExt <- genotypes[(startloc-surrounding):(endloc+surrounding), strain]
        cat(strain, startloc, ":", endloc, posdiff, length(affected), affectedExt)
        if(length(affected) == 3) {
          fill <- NA
          gg <- which(table(affected) == 2)
          if(length(gg) == 1) fill <- names(which(table(affected) == 2))
          fixed <- affected
          fixed[2] <- fill
          newgenotypes[(endloc + startloc) /2, strain] <- fill
          ll <- ll + 1
          cat(" ->", fixed, "")
        }
        if(length(affected) > 3) {
          fill <- NA
          gg <- which(table(affectedExt) == ((2*surrounding) + 2))
          if(length(gg) == 1) fill <- names(which(table(affectedExt) == ((2*surrounding) + 2)))
          newgenotypes[(startloc+1):(endloc-1), strain] <- rep(fill, length((startloc+1):(endloc-1)))
          ll <- ll + 1
          cat(" ", ((2*surrounding) + 2), "->", newgenotypes[(startloc+1):(endloc-1), strain], "")
        }
        cat("\n")
      }
    }
  }
  cat("Fixed", ll, "genotypes\n")
  return(newgenotypes)
}

stats <- calcRecombinations(toNumGeno(genotypes))
genotypes1 <- clean(genotypes, stats$recombinations)

stats1 <- calcRecombinations(toNumGeno(genotypes1))
genotypes2 <- clean(genotypes1, stats1$recombinations)

strains.old <- paste0("BXD", 1:42)        # There are a total of 1848 known recombinations in the 36 older (JAX) BXD set; (~48.1 recombinations per strain)
strains.UTHSC  <- paste0("BXD", 43:103)   # There are a total of 4366 known recombinations in the 53 new (UTHSC) BXD set; (~82.4 recombinations per strain)

mean(stats$recombinationload[strains.old],na.rm=TRUE)
mean(stats1$recombinationload[strains.old],na.rm=TRUE)

mean(stats$recombinationload[strains.UTHSC],na.rm=TRUE)
mean(stats1$recombinationload[strains.UTHSC],na.rm=TRUE)

stats2 <- calcRecombinations(toNumGeno(genotypes2))
genotypes3 <- clean(genotypes2, stats2$recombinations)

mean(stats1$recombinationload[strains.old],na.rm=TRUE)
mean(stats2$recombinationload[strains.old],na.rm=TRUE)

mean(stats1$recombinationload[strains.UTHSC],na.rm=TRUE)
mean(stats2$recombinationload[strains.UTHSC],na.rm=TRUE)


















duplicates <- which(duplicated(numgeno))                                                # Duplicate markers do not add anything

numgeno <- numgeno[-duplicates, ]
physicalmap <- physicalmap[rownames(numgeno),]

cat("Starting with", nrow(numgeno), "of which", length(which(physicalmap[,3] == "Y")), "are selected by rob as markers\n")

chromosomes <- as.character(unique(physicalmap[,"Chr"]))
chrlengths <- unlist(lapply(chromosomes, function(x){ 
  max(physicalmap[which(physicalmap[,"Chr"] == x),"Mb_mm10"])
}))
names(chrlengths) <- chromosomes

plot(x = c(0, length(chromosomes)), y = c(0, max(chrlengths)), t = 'n', xlab = "Chromosome", ylab = "mBp")
lapply(chromosomes){
  
  invisible(NULL)
}


corM <- cor(t(numgeno[physicalmap[,"Chr"] == 1,]), use="pair")
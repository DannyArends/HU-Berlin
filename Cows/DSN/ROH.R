#
# Runs of homozygousity
#

setwd("~/PopStructCow")

# Input files, genotypes 0, 1, 2, phenotypes: breed of individuals
genotypes  = read.table("~/NAS/Cattle/DSN/SNPchip/combined/combined50k_geno_all_uniquesamples_filtered.txt" , sep="\t", header=T, check.names=F)
phenotypes = read.table("~/NAS/Cattle/DSN/SNPchip/combined/combined50k_pheno_all_uniquesamples_filtered.txt", sep="\t", header=T, colClasses="character", quote="", row.names=2)[,-1]

chrs <- gsub("Chr", "", unlist(lapply(strsplit(rownames(genotypes),"_"),"[",1)))
positions <- as.numeric(unlist(lapply(strsplit(rownames(genotypes),"_"),"[",2)))

neworder <- sort(positions, index.return=TRUE)$ix

genotypes <- genotypes[neworder,]
chrs <- chrs[neworder]
positions <- positions[neworder]

orderedgeno <- NULL
for(autosome in 1:28){
  orderedgeno <- rbind(orderedgeno, genotypes[which(chrs == autosome),])
}

# Overwrite the genotypes with the ordered ones
genotypes <- orderedgeno
chrs <- gsub("Chr", "", unlist(lapply(strsplit(rownames(genotypes),"_"),"[",1)))
positions <- as.numeric(unlist(lapply(strsplit(rownames(genotypes),"_"),"[",2)))

# Holstein and DSN
HFind <- rownames(phenotypes)[which(phenotypes[,"breed"] == "Holstein")]
DSNind <- rownames(phenotypes)[which(phenotypes[,"breed"] == "DSN")]

# Is a SNP homozygous (or missing)
isHomozygous <- function(x, homozygous.coding = c(0, 2)){
  if(is.na(x)){ return(TRUE) }
  return(x == homozygous.coding[1] || x == homozygous.coding[2])
}

# Find all runs of homozygousity
findROH <- function(genotype.vector, offset = 0, allow.missing = 0, consecutive.missing = 0, allow.heterozygous = 0, consecutive.heterozygous = 0, verbose=FALSE) {
  index <- 1
  ROHs <- NULL
  while(index < (length(genotype.vector)-1)) {
    if(isHomozygous(genotype.vector[index])) {
      if(verbose) cat("Homozygous genotype at", index, offset, "\n")
      run.index <- index
      n.missing = 0;          # Number of missing values since the start of the run
      con.missing = 0;        # Number of consecutive missing values
      n.heterozygous = 0;     # Number of heterozygous snps since the start of the run
      con.heterozygous = 0;   # Number of consecutive heterozygous snps
      while(run.index < length(genotype.vector)) {
        if (is.na(genotype.vector[run.index])) {
          n.missing <- n.missing + 1        
          con.missing <- con.missing + 1
        } else {
          con.missing <- 0
        }
        if(!isHomozygous(genotype.vector[run.index + 1])){
          n.heterozygous <- n.heterozygous + 1
          con.heterozygous <- con.heterozygous + 1
        }else{
          con.missing <- 0
          con.heterozygous <- 0
        }
        if(n.heterozygous > allow.heterozygous) break;                  # More heterozygous then allowed, break the run
        if(con.heterozygous > consecutive.heterozygous) break;          # More consecutive heterozygous then allowed, break the run
        if(n.missing > allow.missing) break;                            # More missing then allowed, break the run
        if(con.missing > consecutive.missing) break;                    # More consecutive missing then allowed, break the run
        run.index <- run.index + 1
      }
      if(index != run.index) {
        ROHs <- rbind(ROHs, c(index + offset, run.index + offset, n.missing, n.heterozygous, (run.index-index)+1))
      }
    }
    index <- index + 1
  }
  if(length(ROHs) > 0){
    colnames(ROHs) <- c("start", "end", "n.missing", "n.heterozygous", "length")
  }
  return(ROHs)
}

# Greedily get the longest runs above minimum.run.length (no overlap allowed)
longestROHs <- function(ROHs, minimum.run.length = 11, verbose = FALSE) {
  ROHs <- ROHs[sort(ROHs[,"length"], index.return=TRUE, decreasing=TRUE)$ix,]
  ROHs <- ROHs[-which(ROHs[,"length"] < minimum.run.length),]
  index <- 1
  while(index < nrow(ROHs)){
    ROH.start <- ROHs[index,"start"]
    ROH.end <- ROHs[index,"end"]

    idx <- unique(c(which(ROHs[, "end"] >= ROH.start & ROHs[, "end"] <= ROH.end),         #End point is in the region
                    which(ROHs[, "start"] >= ROH.start & ROHs[, "start"] <= ROH.end)))    #Start point is in the region
    idx <- idx[-which(idx == index)]
    if(length(idx) > 0) ROHs <- ROHs[-idx,]
    index <- index + 1
    if(verbose) cat(index, nrow(ROHs), "\n")
  } 
  return(ROHs)
}

# Calculate the ROHs of a single individual
doROH <- function(ind, chrs, genotypes, positions){
  ROHs <- NULL
  offset <- 0
  for(chr in 1:29) {
    snpOnChr <- which(chrs == chr)
    ROHs <- rbind(ROHs, findROH(genotypes[snpOnChr, ind], offset, 5, 3, 3, 2))
    offset <- offset + length(snpOnChr)
  }
  tmp <- longestROHs(ROHs, 20)
  tmp <- cbind(tmp, bplength = NA)
  for(x in 1:nrow(tmp)){
    tmp[x, "bplength"] <- positions[tmp[x,"end"]] - positions[tmp[x,"start"]]
  }
  cat("Done", ind, "\n")
  return(tmp)
}

# Parallel analysis of ROHs in DSN and Holstein
library(parallel)
cl <- makeCluster(getOption("cl.cores", 20))
clusterExport(cl, c("longestROHs", "findROH", "isHomozygous"))
allruns <- parLapply(cl, c(HFind,DSNind), doROH, chrs=chrs, genotypes=genotypes, positions = positions)
stopCluster(cl)
names(allruns) <- c(HFind,DSNind)
save(allruns, file="~/PopStructCow/AllROH.Rdata")

#Get the ROH as a genotype matrix (0, 1)
getROHmatrix <- function(allruns, genotypes, individuals, size=4000000){
  allruns.subset <- allruns[individuals]
  ROH.geno <- matrix(0, nrow(genotypes), length(individuals), dimnames=list(rownames(genotypes), individuals))
  for(ind in c(individuals)) {
    ROHs <- allruns[[ind]][which(allruns[[ind]][,"bplength"] > size),]
    if(length(ROHs) == 0) next; # Nothing
    if(length(ROHs) == 6){ # Vector
      ROH.geno[ROHs["start"]:ROHs["end"], ind] <- ROH.geno[ROHs["start"]:ROHs["end"], ind] + 1
    }else{ # Matrix
      for(x in 1:nrow(ROHs)){
        ROH.geno[ROHs[x,"start"]:ROHs[x,"end"], ind] <- ROH.geno[ROHs[x,"start"]:ROHs[x,"end"], ind] + 1
      }
    }
  }
  return(ROH.geno)
}

## Sumarized positions for in the plot
sumpos <- positions
pos.add <- 0
for(chr in unique(chrs)){
  chr.length <- max(positions[which(chrs==chr)],na.rm=TRUE)
  sumpos[which(chrs==chr)] <- sumpos[which(chrs==chr)] + pos.add
  pos.add <- pos.add + chr.length + 10000000
}
sumpos <- sumpos / 1000000

### minimum 4 MB ROH

ROH.HF.4mb.geno <- getROHmatrix(allruns, genotypes, HFind, 4000000)
ROH.DSN.4mb.geno <- getROHmatrix(allruns, genotypes, DSNind, 4000000)

ROH.HF.4mb.vector <- apply(ROH.HF.4mb.geno,1,sum)/ ncol(ROH.HF.4mb.geno)
cutoff.HF.4mb <- median(ROH.HF.4mb.vector) + 2* sd(ROH.HF.4mb.vector)

ROH.DSN.4mb.vector <- apply(ROH.DSN.4mb.geno,1,sum)/ ncol(ROH.DSN.4mb.geno)
cutoff.DSN.4mb <- median(ROH.DSN.4mb.vector) + 2* sd(ROH.DSN.4mb.vector)

### minimum 8 MB ROH

ROH.HF.8mb.geno <- getROHmatrix(allruns, genotypes, HFind, 8000000)
ROH.DSN.8mb.geno <- getROHmatrix(allruns, genotypes, DSNind, 8000000)

ROH.HF.8mb.vector <- apply(ROH.HF.8mb.geno,1,sum)/ ncol(ROH.HF.8mb.geno)
cutoff.HF.8mb <- median(ROH.HF.8mb.vector) + 2* sd(ROH.HF.8mb.vector)

ROH.DSN.8mb.vector <- apply(ROH.DSN.8mb.geno,1,sum)/ ncol(ROH.DSN.8mb.geno)
cutoff.DSN.8mb <- median(ROH.DSN.8mb.vector) + 2* sd(ROH.DSN.8mb.vector)

### minimum 16 MB ROH

ROH.HF.16mb.geno <- getROHmatrix(allruns, genotypes, HFind, 16000000)
ROH.DSN.16mb.geno <- getROHmatrix(allruns, genotypes, DSNind, 16000000)

ROH.HF.16mb.vector <- apply(ROH.HF.16mb.geno,1,sum)/ ncol(ROH.HF.16mb.geno)
cutoff.HF.16mb <- median(ROH.HF.16mb.vector) + 2* sd(ROH.HF.16mb.vector)

ROH.DSN.16mb.vector <- apply(ROH.DSN.16mb.geno,1,sum)/ ncol(ROH.DSN.16mb.geno)
cutoff.DSN.16mb <- median(ROH.DSN.16mb.vector) + 2* sd(ROH.DSN.16mb.vector)

### ROH plots
png("ROH_genomewide.png", width=1024, height=600)

op <- par(mfrow=c(3,1))
colz <- c("blue", "orange")[as.numeric(as.numeric(chrs) %% 2 == 0)+1]

maxy = max(ROH.HF.4mb.vector, ROH.DSN.4mb.vector) * 1.15

plot(x=c(0, max(sumpos,na.rm=TRUE)), y=c(-maxy, maxy), t = 'n', main="ROH (4mb)", ylab="Poportion in ROH", xlab="Position")
points(x=sumpos, y=ROH.HF.4mb.vector, col=colz, pch=19, cex=0.4)
abline(h = cutoff.HF.4mb)
points(x=sumpos, y=-ROH.DSN.4mb.vector, col=colz, pch=19, cex=0.4)
abline(h = -cutoff.DSN.4mb)

maxy = max(ROH.HF.8mb.vector, ROH.DSN.8mb.vector) * 1.15

plot(x=c(0, max(sumpos,na.rm=TRUE)), y=c(-maxy, maxy), t = 'n', main="ROH (8mb)", ylab="Poportion in ROH", xlab="Position")
points(x=sumpos, y=ROH.HF.8mb.vector, col=colz, pch=19, cex=0.4)
abline(h = cutoff.HF.8mb)
points(x=sumpos, y=-ROH.DSN.8mb.vector, col=colz, pch=19, cex=0.4)
abline(h = -cutoff.DSN.8mb)

maxy = max(ROH.HF.16mb.vector, ROH.DSN.16mb.vector) * 1.15

plot(x=c(0, max(sumpos,na.rm=TRUE)), y=c(-maxy, maxy), t = 'n', main="ROH (16mb)", ylab="Poportion in ROH", xlab="Position")
points(x=sumpos, y=ROH.HF.16mb.vector, col=colz, pch=19, cex=0.4)
abline(h = cutoff.HF.16mb)
points(x=sumpos, y=-ROH.DSN.16mb.vector, col=colz, pch=19, cex=0.4)
abline(h = -cutoff.DSN.16mb)

dev.off()
### Total length of genome covered by SNPs

chr.length.total <- 0
for(chr in unique(chrs)){
  chr.length.total <- chr.length.total + (max(positions[which(chrs==chr)],na.rm=TRUE) - min(positions[which(chrs==chr)],na.rm=TRUE))
}

### Sum of length covered by ROH in an individual
getROHlength <- function(allruns, ind, size=4000000) {
  ROHs <- allruns[[ind]][which(allruns[[ind]][,"bplength"] > size),]
  if(length(ROHs) == 0) return(0)
  if(length(ROHs) == 6) return(ROHs["bplength"])
  return(sum(ROHs[,"bplength"]))
}

Frohs.4mb.HF <- Frohs.8mb.HF <- Frohs.16mb.HF <- NULL
for(ind in HFind){
  Frohs.4mb.HF <- c(Frohs.4mb.HF, getROHlength(allruns, ind, 4000000) / chr.length.total)
  Frohs.8mb.HF <- c(Frohs.8mb.HF, getROHlength(allruns, ind, 8000000) / chr.length.total)
  Frohs.16mb.HF <- c(Frohs.16mb.HF, getROHlength(allruns, ind, 16000000) / chr.length.total)
}
HF.4mb.hist <- hist(Frohs.4mb.HF, breaks=seq(0,1,0.01), plot=FALSE)
HF.8mb.hist <- hist(Frohs.8mb.HF, breaks=seq(0,1,0.01), plot=FALSE)
HF.16mb.hist <- hist(Frohs.16mb.HF, breaks=seq(0,1,0.01), plot=FALSE)

Frohs.4mb.DSN <- Frohs.8mb.DSN <- Frohs.16mb.DSN <- NULL
for(ind in DSNind){
  Frohs.4mb.DSN <- c(Frohs.4mb.DSN, getROHlength(allruns, ind, 4000000) / chr.length.total)
  Frohs.8mb.DSN <- c(Frohs.8mb.DSN, getROHlength(allruns, ind, 8000000) / chr.length.total)
  Frohs.16mb.DSN <- c(Frohs.16mb.DSN, getROHlength(allruns, ind, 16000000) / chr.length.total)
}
DSN.4mb.hist <- hist(Frohs.4mb.DSN, breaks=seq(0,1,0.01), plot=FALSE)
DSN.8mb.hist <- hist(Frohs.8mb.DSN, breaks=seq(0,1,0.01), plot=FALSE)
DSN.16mb.hist <- hist(Frohs.16mb.DSN, breaks=seq(0,1,0.01), plot=FALSE)

### Plot the F(ROH)
png("F_ROH.png", width=800, height=800)

xpos <- HF.4mb.hist$breaks[-length(HF.4mb.hist$breaks)] + 0.5 * diff(HF.4mb.hist$breaks)

plot(c(0,0.3), c(0,0.5), t = 'n', ylab="Poportion in ROH", xlab="Inbreeding coefficient F(ROH)")
points(x = xpos, y = HF.4mb.hist$counts / sum(HF.4mb.hist$counts), pch=19, t='l')
points(x = xpos, y = HF.4mb.hist$counts / sum(HF.4mb.hist$counts), pch=19, t='p')
points(x = xpos, y = HF.8mb.hist$counts / sum(HF.8mb.hist$counts), pch=17, t='l')
points(x = xpos, y = HF.8mb.hist$counts / sum(HF.8mb.hist$counts), pch=17, t='p')

points(x = xpos , y = DSN.4mb.hist$counts / sum(DSN.4mb.hist$counts), col=2, pch=19, t='l')
points(x = xpos , y = DSN.4mb.hist$counts / sum(DSN.4mb.hist$counts), col=2, pch=19, t='p')
points(x = xpos , y = DSN.8mb.hist$counts / sum(DSN.8mb.hist$counts), col=2, pch=17, t='l')
points(x = xpos , y = DSN.8mb.hist$counts / sum(DSN.8mb.hist$counts), col=2, pch=17, t='p')

legend("topright", c("Holstein roh>4mb", "Holstein roh>8mb","DSN roh>4mb","DSN roh>8mb"), pch=c(19,17,19,17), col=c(1,1,2,2))

dev.off()

### Plot the ROH length versus the number of ROH segments
### Sum of length covered by ROH in an individual
getROHs <- function(allruns, ind, size=4000000) {
  ROHs <- allruns[[ind]][which(allruns[[ind]][,"bplength"] > size),]
  return(ROHs)
}

ROH.segments.DSN <- NULL
for(ind in DSNind){
  ROH <- getROHs(allruns, ind)
  if(length(ROH) == 0){
    ROH.segments.DSN <- rbind(ROH.segments.DSN, c(0,0))
  }else if(length(ROH) == 6){
    ROH.segments.DSN <- rbind(ROH.segments.DSN, c(1,ROH["bplength"]))
  }else{
    ROH.segments.DSN <- rbind(ROH.segments.DSN, c(nrow(ROH), sum(ROH[,"bplength"])))
  }
}

ROH.segments.HF <- NULL
for(ind in HFind){
  ROH <- getROHs(allruns, ind)
  if(length(ROH) == 0){
    ROH.segments.HF <- rbind(ROH.segments.HF, c(0,0))
  }else if(length(ROH) == 6){
    ROH.segments.HF <- rbind(ROH.segments.HF, c(1,ROH["bplength"]))
  }else{
    ROH.segments.HF <- rbind(ROH.segments.HF, c(nrow(ROH), sum(ROH[,"bplength"])))
  }
}

ROH.segments.DSN[,2]  <- ROH.segments.DSN[,2] / 1000000
ROH.segments.HF[,2]  <- ROH.segments.HF[,2] / 1000000

xmax <- max(c(ROH.segments.HF[,2], ROH.segments.DSN[,2]))
ymax <- max(c(ROH.segments.HF[,1], ROH.segments.DSN[,1]))

png("ROH_segments.png", width=800, height=800)

plot(c(0,xmax), c(0,ymax), t ='n', xlab="Covered genome length (Mb)", ylab="Number of ROH segements")
points(ROH.segments.HF[,2], ROH.segments.HF[,1], col=1, pch=19,cex=0.7)
points(ROH.segments.DSN[,2], ROH.segments.DSN[,1], col=2, pch=19,cex=0.7)

legend("bottomright", c("Holstein roh>4mb","DSN roh>4mb"), pch=c(19,19), col=c(1,2))

dev.off()


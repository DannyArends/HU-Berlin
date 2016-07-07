#
# Server script to extract SNPs for a region
#

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); } #q("no") }
}

samples <- c("ext_L7254", "ext_L7256", "ext_L7257", "ext_L7255", "ext_L7258")
names(samples) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")

# Start of analysis
setwd("~/DIFE/")
regions <- read.table("regions_bfmi.txt", sep = "\t",header=TRUE)   # Start loading the BFMi regions

# Use tabix to subset the big vcf.gz file
for(i in 1:nrow(regions)){
  outfile <- paste0("~/DIFE/analysis/", regions[i,"id"], ".txt")
  cmd <- paste0("tabix -h /home/arends/NAS/Mouse/DNA/Sequencing/DifeMouse/RAW/ALL_variants.vcf.gz chr", regions[i,"chr"], ":", regions[i,"start"], "-", regions[i,"end"], " > ", outfile)
  execute(cmd, outfile)
}

# Plot all regions for a certain chromosome
plotRegions <- function(regions, chr = 1){
  chrsubset <- regions[which(regions[,"chr"] == chr),]
  mY <- max(chrsubset[, "end"], na.rm=TRUE)
  mX <- nrow(chrsubset)
  plot(c(0,mX), c(0,mY), t = 'n', yaxt='n')
  for(i in 1:nrow(chrsubset)){
    points(c(i,i), c(chrsubset[i,"start"], chrsubset[i, "end"]), t = 'l', col=i,lwd=2)
  }
  axis(2, seq(0, mY, 10000000), seq(0, mY, 10000000)/ 1000000, las=2)
}

plotRegions(regions, 3)

loadRegion <- function(name = "HUregion1", loc = "~/DIFE/analysis/") {
  header <- strsplit(readLines(paste0(loc, name, ".txt"), n = 133)[133], "\t")[[1]]
  header[1] <- "CHROM"
  regiondata <- read.table(paste0(loc, name, ".txt"),sep="\t")
  colnames(regiondata) <- header
  idx <- which(colnames(regiondata) %in% samples)
  colnames(regiondata)[idx] <- names(samples)[match(colnames(regiondata)[idx], samples)]

  snpnames <- as.character(regiondata[,"ID"])
  noID <- which(snpnames == ".")
  snpnames[noID] <- paste0("N-", gsub("chr", "", regiondata[noID, "CHROM"]), ":", regiondata[noID, "POS"])
  rownames(regiondata) <- snpnames
  return(regiondata)
}

onlyGeno <- function(regiondata){
  geno <- t(apply(regiondata[,names(samples)], 1, function(x){
    unlist(lapply(strsplit(x, ":"), "[", 1))
  }))
  rownames(geno) <- rownames(regiondata)
  return(geno)
}

filterGeno <- function(genotypes, mfilter, naIsOK = FALSE){
  dontcare <- names(which(is.na(mfilter)))
  geno <- genotypes[, -which(colnames(genotypes) %in% dontcare)]
  mfilter <- mfilter[-which(is.na(mfilter))]
  
  matches <- NULL
  for(x in 1:nrow(geno)){
    if(all(geno[x, names(mfilter)] == mfilter)){
      matches <- cbind(matches, x)
    }
  }
  return(genotypes[matches,])
}

regiondata <- loadRegion()
genotypes <- onlyGeno(regiondata)

mfilterS1 <- c(NA, "1/1", "0/0", NA, NA); names(mfilterS1) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")
mfilterS2 <- c(NA, "0/0", "1/1", NA, NA); names(mfilterS2) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")

S1S2diff <- rbind(filterGeno(genotypes, mfilterS1), filterGeno(genotypes, mfilterS2))


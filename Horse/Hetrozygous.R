# Analysis of Kabadiner horse SNP and performance data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015
setwd("E:/Horse/DNA/Kabadiner/")
                                                                                                  # Read in the data from the Mega Muga
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

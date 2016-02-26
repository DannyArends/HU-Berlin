
setwd("E:/Mouse/DNA/MegaMuga/")

# Load the Allele Transmission Biased regions
regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)

allregions <- rbind(regionsMat, regionsPat)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross/Analysis")
combinedASE <- read.table("combined_ASE_3tissues.txt", sep ="\t", header=TRUE)

for(r in 1:nrow(allregions)) {
  chr <- allregions[r, "Chr"]
  rstart <- allregions[r, "Start"]
  rend <- allregions[r, "Stop"]
  ii <- which(combinedASE[,"CHROM"] == chr & combinedASE[,"POS"] > rstart & combinedASE[,"POS"] < rend)
  cat(r," ", length(ii), "\n")
  if(r == 35) break;
}

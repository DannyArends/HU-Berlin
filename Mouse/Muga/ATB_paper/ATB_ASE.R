
setwd("E:/Mouse/DNA/MegaMuga/")

# Load the Allele Transmission Biased regions
regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)

allregions <- rbind(regionsMat, regionsPat)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross/Analysis")
combinedASE <- read.table("combined_ASE_3tissues.txt", sep ="\t", header=TRUE)

classes <- c("diff", "mB6N", "mBFMI")

significant <- vector("list", 3)
cnt <- 1
for(x in classes){
  significant[[cnt]] <- c(significant[[cnt]], apply(read.csv(paste0("snp_stats_Gonadal Fat_ASE_",x,".txt"), sep ="\t", colClasses="character")[,c("CHROM", "POS", "REF", "ALT", "transcript_id")], 1, paste0, collapse="_"))
  significant[[cnt]] <- c(significant[[cnt]], apply(read.csv(paste0("snp_stats_Liver_ASE_",x,".txt"), sep ="\t", colClasses="character")[,c("CHROM", "POS", "REF", "ALT", "transcript_id")], 1, paste0, collapse="_"))
  significant[[cnt]] <- c(significant[[cnt]], apply(read.csv(paste0("snp_stats_Quadriceps_ASE_",x,".txt"), sep ="\t", colClasses="character")[,c("CHROM", "POS", "REF", "ALT", "transcript_id")], 1, paste0, collapse="_"))
  cnt <- cnt +1
}

Diff <- combinedASE[which(combinedASE[,"ID"] %in% unique(significant[[1]])),]
cat("Different", nrow(Diff), "\n")
mB6N <- combinedASE[which(combinedASE[,"ID"] %in% unique(significant[[2]])),]
cat("Different", nrow(mB6N), "\n")
mBFMI <- combinedASE[which(combinedASE[,"ID"] %in% unique(significant[[3]])),]
cat("Different", nrow(mBFMI), "\n")

matBFMI <- c(c("BFMI860.12xB6N_4425_GonadalFat", "BFMI860.12xB6N_4426_GonadalFat", "BFMI860.12xB6N_4427_GonadalFat"), 
c("BFMI860.12xB6N_5073_Liver","BFMI860.12xB6N_5074_Liver","BFMI860.12xB6N_5075_Liver"), 
c("BFMI860.12xB6N_5080_Quadriceps","BFMI860.12xB6N_5081_Quadriceps","BFMI860.12xB6N_5079_Quadriceps"))

matB6N <- c(c("B6NxBFMI860.12_4424_GonadalFat","B6NxBFMI860.12_4422_GonadalFat","B6NxBFMI860.12_4423_GonadalFat"), 
c("B6NxBFMI860.12_5072_Liver","B6NxBFMI860.12_5070_Liver","B6NxBFMI860.12_5071_Liver"), 
c("B6NxBFMI860.12_5076_Quadriceps","B6NxBFMI860.12_5078_Quadriceps","B6NxBFMI860.12_5077_Quadriceps"))

# THIS DOES NOT WORK, we cannot use the ASE data since we looed for differential ASE between matBFMI and matB6N
# Additionally, we have 3 individuals sumamrized into 1 number per tissue which is wrong (54, 46 44) -> Average below 50, but not consistent, etc
for(r in 1:nrow(allregions)) {
  chr <- allregions[r, "Chr"]
  rstart <- allregions[r, "Start"]
  rend <- allregions[r, "Stop"]
  Diffii <- which(Diff[,"CHROM"] == chr & Diff[,"POS"] > rstart & Diff[,"POS"] < rend)
  mB6Nii <- which(mB6N[,"CHROM"] == chr & mB6N[,"POS"] > rstart & mB6N[,"POS"] < rend)
  mBFMIii <- which(mBFMI[,"CHROM"] == chr & mBFMI[,"POS"] > rstart & mBFMI[,"POS"] < rend)
  cat(r," ", length(Diffii)," ", length(mB6Nii)," ", length(mBFMIii), "\n")
  
  Diffii <- unique(c(Diffii, mB6Nii, mBFMIii))
  
  if(length(Diffii) > 0){
    
    Diffsubset <- cbind(Diff[Diffii, 1:11], 
        t(apply(Diff[Diffii, matBFMI],1, function(x){ return(c(mean(x[1:3]),mean(x[4:6]),mean(x[7:9]))) })),
        t(apply(Diff[Diffii, matB6N],1, function(x){ return(c(mean(x[1:3]),mean(x[4:6]),mean(x[7:9]))) })),
        BFMI = unique(na.omit(Diff[Diffii, "BFMI860.12_6344_call_GonadalFat"], Diff[Diffii, "BFMI860.12_5067_call_Liver"], Diff[Diffii, "BFMI860.12_6340_call_Quadriceps"]))
    )
    colnames(Diffsubset)[12:17] <- c("matBFMI","matBFMI","matBFMI","matB6N","matB6N","matB6N")
    
    for(x in 1:nrow(Diffsubset)){
      if(Diffsubset[x,"BFMI"] == "Alt"){
        nii <- which(Diffsubset[x, 12:17] < 50)
        mii <- which(Diffsubset[x, 12:17] > 50)
        cat(length(nii), length(mii), "\n")
      }
    }
  }
    
  if(r == 60) break
}



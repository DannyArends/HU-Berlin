setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
results <- read.table("ChiSqScores_Incompatible.txt", sep="\t", check.names=FALSE)
rr <- apply(results,1,as.numeric)
rownames(rr) <- colnames(rr)
results <- rr

# Load the Allele Transmission Biased regions
regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)

allRegions <- rbind(cbind(regionsMat, origin = "MAT"), cbind(regionsPat, origin = "PAT"))
allRegions <- allRegions[with(allRegions, order(Chr, Start)), ]

regionnames <- paste0(allRegions[,"Chr"],":", 
          round(allRegions[,"Start"] / 1000000, 0),"-",
          round(allRegions[,"Stop"] / 1000000, 0), " ", 
          allRegions[,"origin"])

rownames(allRegions) <- regionnames

allRegions["8:62-63 PAT", "Prefered.Allele"] <- "B6N"
allRegions["9:88-91 MAT", "Prefered.Allele"] <- "B6N"
allRegions["11:13-13 MAT", "Prefered.Allele"] <- "BFMI"


# Generate an overview plot
nTests <- (ncol(results) * ncol(results)) / 2
LODscores <- -log10(pchisq(results, 1, lower.tail=FALSE))
threshold <- -log10(0.05 / nTests)

significant <-  names(which(apply(LODscores,1,max,na.rm=TRUE) > threshold))

allRegions <- allRegions[significant, ]
results <- results[significant,significant]



bfmipref <- which(allRegions[,"Prefered.Allele"] == "BFMI")
b6npref <- which(allRegions[,"Prefered.Allele"] == "B6N")

pat_bfmi_snp <- read.table(file="Analysis/PAT_BFMI_NSyn_Filtered.txt", sep = "\t")
mat_bfmi_snp <- read.table(file="Analysis/MAT_BFMI_NSyn_Filtered.txt", sep = "\t")
pat_b6_snp <- read.table(file="Analysis/PAT_B6N_NSyn_Filtered.txt", sep = "\t")
mat_b6_snp <- read.table(file="Analysis/MAT_B6N_NSyn_Filtered.txt", sep = "\t")

for(x in rownames(results)){
  genes <- NULL
  cat("------------------------\n")
  fp <- paste0("Analysis/", allRegions[x, "origin"], "_", allRegions[x, "Prefered.Allele"], "_NSyn_Filtered.txt")
  geneTableO <- read.table(file = fp, sep = "\t", header=TRUE)
  splitted <- strsplit(x, ":")[[1]]
  chrO <- splitted[1]
  posO <- strsplit(splitted[2], " ")[[1]][1]
  posOs <- (as.numeric(strsplit(posO, "-")[[1]][1])-1) * 1000000
  posOe <- (as.numeric(strsplit(posO, "-")[[1]][2])+1) * 1000000
  inR <- which(geneTableO[, "chromosome_name"] == chrO & geneTableO[, "start_position"] > posOs & geneTableO[, "end_position"] < posOe)
  if(length(inR) > 0) {
    genes <- as.character(geneTableO[inR, "mgi_symbol"])
    for(g in genes){
      cat(x, g, " ", "*\n")
    }
    for(y in names(which(results[x,] > threshold))){
      fp2 <- paste0("Analysis/", allRegions[y, "origin"], "_", allRegions[y, "Prefered.Allele"], "_NSyn_Filtered.txt")
      geneTableR <- read.table(file = fp, sep = "\t", header=TRUE)
      splitted <- strsplit(y, ":")[[1]]
      chrR <- splitted[1]
      posR <- strsplit(splitted[2], " ")[[1]][1]
      posRs <- (as.numeric(strsplit(posR, "-")[[1]][1])-1) * 1000000
      posRe <- (as.numeric(strsplit(posR, "-")[[1]][2])+1) * 1000000
      inR2 <- which(geneTableO[, "chromosome_name"] == chrR & geneTableO[, "start_position"] > posRs & geneTableO[, "end_position"] < posRe)
      if(length(inR2) > 0) {
        genes2 <- as.character(geneTableR[inR2, "mgi_symbol"])
        for(g in genes2){
          cat(y, g, " ", "+\n")
        }
      }
    }
  }
  cat("------------------------\n")
}


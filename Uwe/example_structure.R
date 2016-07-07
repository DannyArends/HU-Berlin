
setwd("E:/UWE")

# Raw SNP data
snpdata <- read.csv("Ziegen_HU-Berlin_Matrix.txt", sep="\t", skip = 9, header=TRUE, check.names=FALSE, row.names=1, na.strings=c("NA", "", "--"))
# SNP information (Allele, locations etc)
snpinfo <- read.table("snpinfo.txt", sep="\t", na.strings=c("", "NA", "N.D."))

snpdata <- snpdata[rownames(snpinfo), ]                         # Take only the SNPs for which we have information (2000 here)
snpdata <- snpdata[,-which(colnames(snpdata) == "DN 2")]        # Throw away the duplicate individual because it confuses STRUCTURE

### Split the SNP alleles from C/T to [C, T]
snpAlleles <- lapply(strsplit(as.character(snpinfo[,"allele"]), ""), "[", c(1,3))

### Transform every row of the SNP matrix from AC, CC to 1 and 2
numsnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames = list(rownames(snpdata), colnames(snpdata)))
for(x in 1:length(snpAlleles)) {
  if(!is.na(snpinfo[x, "reference"]) && snpAlleles[[x]][1] !=  snpinfo[x, "reference"]){  # C/T while reference is T, so flip it around
    snpAlleles[[x]] <- snpAlleles[[x]][2:1]
  }
  g1 <- paste(snpAlleles[[x]][1], snpAlleles[[x]][1],sep="")
  g2a <- paste(snpAlleles[[x]][1], snpAlleles[[x]][2],sep="")
  g2b <- paste(snpAlleles[[x]][2], snpAlleles[[x]][1],sep="")
  g3 <- paste(snpAlleles[[x]][2], snpAlleles[[x]][2],sep="")
  if(!all(as.character(unlist(snpdata[x,])) %in% c(g1,g2a,g2b,g3, NA))) stop("Nope")
  numsnpdata[x, which(snpdata[x, ] == g1)] <- 1     # A
  numsnpdata[x, which(snpdata[x, ] == g2a)] <- 2    # H
  numsnpdata[x, which(snpdata[x, ] == g2b)] <- 2    # H
  numsnpdata[x, which(snpdata[x, ] == g3)] <- 3     # B
}
write.table(numsnpdata, "snps_numeric_NO_DN2.txt", sep="\t", quote=FALSE)            # Write the numeric files

### Write out the data for STRUCTURE
numsnpdata <- t(numsnpdata)
structGeno <- NULL
for(x in 1:nrow(numsnpdata)){
  gg <- rbind(rep(NA, ncol(numsnpdata)), rep(NA, ncol(numsnpdata)))
  a1 <- which(numsnpdata[x,] == 1)
  a2 <- which(numsnpdata[x,] == 2)
  a3 <- which(numsnpdata[x,] == 3)
  gg[1, a1] <- 0; gg[2, a1] <- 0                                                      # Pretty inefficient, but it will do the trick
  gg[1, a2] <- 0; gg[2, a2] <- 1
  gg[1, a3] <- 1; gg[2, a3] <- 1
  gg[is.na(gg)] <- 9
  structGeno <- rbind(structGeno, gg)
}

rownames(structGeno) <- gsub(" ","", unlist(lapply(rownames(numsnpdata), rep, 2)))    # Spaces are not allowed in sample names
colnames(structGeno) <- colnames(numsnpdata)
write.table(structGeno, file="snps_structure_NO_DN2.txt", sep = "\t")                 # Save the genotypes to disk

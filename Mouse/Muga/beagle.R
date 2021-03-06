# Phasing the MUGA SNP array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sep, 2014
# first written Sep, 2014

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                              # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

setwd("E:/Mouse/DNA/MegaMuga/")                                                             # Read in the data from the Mega Muga
load("MM_snps.Rdata")                                                                       # JAX data
MM_snps[,"SNP_ID"] <- gsub(".", "-", MM_snps[,"SNP_ID"], fixed = TRUE)                      # Fix the . to - error in the MM_snps from JAX

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))                  # Genetic map (created from JAX data)
genotypes <- read.table("Analysis/genotypesRAW.txt", sep="\t", check.names=FALSE)           # Genotypes measured on the MUGA array

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })        # Missing amount of genotype data per individual
genotypes <- genotypes[,-which(missingPerInd==100)]                                         # Remove individuals which have NO genotypes

alleles <- NULL                                                                             # Create the alleles
for(x in 1:nrow(MM_snps)){
  instr <- regexpr("\\[.+\\]", MM_snps[x,"Sequence"],perl=TRUE)
  values <- unlist(strsplit(substr(MM_snps[x,"Sequence"], (instr[1]+1), (instr[1]+attr(instr,"match.length")-2)),"/"))
  alleles <- rbind(alleles, values)
}

rownames(alleles) <- MM_snps[,"SNP_ID"]
alleles <- alleles[rownames(genotypes),]                                                    # Sort the alleles so they match

phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),]       # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]
P  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]

hetroPerc <- NULL
for(ind in P){ hetroPerc <- c(hetroPerc, sum(genotypes[, ind] == "H",na.rm=TRUE)/nrow(genotypes) * 100) }

toSNPgenotype <- function(allele, alleles, sep = " "){                                      # SNP genotype function, from A, H, B to the actual SNPs
  if(is.na(allele)) return(paste("?", "?", sep=sep))
  if(allele == "A"){ return(paste(alleles[1], alleles[1], sep=sep)); }
  if(allele == "H"){ return(paste(alleles[1], alleles[2], sep=sep)); }
  if(allele == "B"){ return(paste(alleles[2], alleles[2], sep=sep)); }
}

checkAlleles <- function(pA, mA, cA){                                                       # Check alleles for segregation errors
  if(is.na(cA)) return(TRUE)
  if(is.na(pA) && is.na(mA)) return(TRUE)
  if(is.na(pA) && mA == "A") return(cA != "B")
  if(is.na(pA) && mA == "H") return(TRUE)
  if(is.na(pA) && mA == "B") return(cA != "A")
  if(is.na(mA) && pA == "A") return(cA != "B")
  if(is.na(mA) && pA == "H") return(TRUE)
  if(is.na(mA) && pA == "B") return(cA != "A")
  if(pA != "H" && pA == mA) return(cA==pA)
  if(pA == "H" && mA == "A") return(cA != "B")
  if(pA == "H" && mA == "B") return(cA != "A")
  if(mA == "H" && pA == "A") return(cA != "B")
  if(mA == "H" && pA == "B") return(cA != "A")
  if(pA == "A" && mA == "B") return(cA=="H")
  if(pA == "B" && mA == "A") return(cA=="H")
  return(TRUE)
}

errors <- NULL
for(ind in c(F1, F2)){                                                                         # Phase all the children
  s <- proc.time()
  err <- 0
  vID <- which(colnames(genotypes) %in% phenotypes[ind,"Vater"])
  mID <- which(colnames(genotypes) %in% phenotypes[ind,"Mutter"])

  datamatrix <- NULL
  if(length(vID) > 0 && length(mID) > 0){                                                     # We have both parents in the dataset
    for(x in 1:nrow(genotypes)){
      if(checkAlleles(genotypes[x, vID], genotypes[x, mID], genotypes[x, ind])){
        datamatrix <- rbind(datamatrix, c("M",rownames(genotypes)[x],toSNPgenotype(genotypes[x, vID], alleles[x,]), toSNPgenotype(genotypes[x, mID], alleles[x,]), toSNPgenotype(genotypes[x, ind], alleles[x,])))
      }else{
        err <- err + 1
        #cat("genotype error",ind,"at:", rownames(genotypes)[x],"->", as.character(genotypes[x, vID]), as.character(genotypes[x, mID]), as.character(genotypes[x, ind]), "\n")
      }
    }
  }else{
    cat("Not all parents are present in the data for", ind, "\n")
  }

  input  <- paste0("Analysis/phase/", ind, ".bgl")
  output <- paste0("Analysis/phase/", ind)
  errors <- c(errors, err)

  write.table(datamatrix, file=input, sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)
  system(paste0("java -Xmx1000m -jar beagle.jar redundant=TRUE trios=", input, " missing=? out=", output))
  cat("--- Done", ind, "found", err, "errors in", (proc.time() -s)[3], "second\n")
}
png(file="Analysis/OffspringGenotypeMismatches.png")
  plot(errors, pch=19)
dev.off()

c(F1, F2)[which(errors > 500)]                                                              # individuals that have 500+ mismatches between parents and children
### [1] "33310233" "6661965"  "6662155"  "6662156" 

phasedGeno <- matrix("", nrow(genotypes), ncol(genotypes))                                  # Empty matrix for the phased genotypes
rownames(phasedGeno) <- rownames(genotypes) ; colnames(phasedGeno) <- colnames(genotypes)   # Set names

datamatrix <- NULL
for(m in 1:nrow(genotypes)){                                                                # Phase the parent data (we cannot use trio information)
  mrow <- c("M", rownames(genotypes)[m])
  for(ind in P){ 
    mrow <- c(mrow, toSNPgenotype(genotypes[m, ind], alleles[m,],""))
  } 
  datamatrix <- rbind(datamatrix, mrow)
}
write.table(datamatrix, file="Analysis/phase/parents.bgl", sep=" ",quote=FALSE, col.names=FALSE, row.names=FALSE)
system(paste0("java -Xmx1000m -jar beagle.jar redundant=TRUE unphased=Analysis/phase/parents.bgl missing=? out=Analysis/phase/parents"))

phaseddata <- read.table(gzfile(paste0("Analysis/phase/parents.parents.bgl.phased.gz")), header=TRUE, colClasses="character", row.names=2)
colnames(phaseddata) <- c("M", P)

for(ind in P){ phasedGeno[rownames(phaseddata), ind] <- phaseddata[, ind] }

# Get the F1 and F2 data into the matrix
for(ind in c(F1,F2)){
  phaseddata <- read.table(gzfile(paste0("Analysis/phase/", ind,".",ind,".bgl.phased.gz")), header=TRUE, colClasses="character", row.names=2)
  phasedGeno[rownames(phaseddata), ind] <-  apply(phaseddata[,6:7],1,function(x){paste0(x[1],x[2]); })
}

write.table(phasedGeno, "Analysis/genotypesPhasedBeagle.txt", sep="\t")

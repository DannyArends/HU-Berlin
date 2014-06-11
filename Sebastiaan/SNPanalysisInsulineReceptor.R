# SNPanalysisInsulineReceptor.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Analysis of data from Sebastiaan

setwd("D://Sebastiaan_mouse")

chromosomes <- c(1:19, "M", "X", "Y")
rowheader   <- 1:6
combinedFN  <- "S-1versusBFMI.txt"

aa <- lapply(chromosomes, function(x){
  chrdata <- read.csv(paste0("Mm_chr", x, "_snp_genotypes_cleaned.txt"), colClasses = "character", sep="\t", skip = 4)

  chr <- chrdata[c(6:nrow(chrdata)),]
  colnames(chr)[1:6] <- as.character(unlist(chrdata[5,1:6]))
#  colnames(chr)[7:ncol(chr)] <- as.character(unlist(chrdata[4,7:ncol(chr)]))

  inconsistent <- which(chr[,"BFMI861.S1"] != chr[,"BFMI861.S1.1"])                         # The inconsistent SNPs _S1 and -S1
  inconsistent <- c(inconsistent, which(chr[,"BFMI861.S2"] != chr[,"BFMI861.S2.1"]))        # The inconsistent SNPs _S2 and -S2
  inconsistent <- c(inconsistent, which(chr[,"BFMI861.S1"] == "00"))                       # The SNPs not assessed in BFMI861.S1
  inconsistent <- unique(inconsistent)                                                      # Get only the unique ones
  
  chr <- chr[-inconsistent, ]                                                                                 # Throw away the nonconsistent SNPs
  
  nonDiabetic <- c("BFMI852", "BFMI856", "BFMI860.12", "BFMI860.12.1" , "BFMI860.S2", "BFMI861.S2.1")         # nonDiabetic BFMI mice
  isDiabetic <- c("BFMI861.S1")                                                             # Diabetic strain

  if(nrow(chr) > 0){
    dSNP <- NULL
    for(snp in 1:nrow(chr)){                                                                # For all the snp
      if(!(chr[snp, isDiabetic] %in% chr[snp, nonDiabetic])){ dSNP <- c(dSNP, snp) }        # If the diabetic SNP is not in the non Diabetic group, Add the SNP
    }

    snpOUT <- cbind(chr[dSNP, rowheader], chr[dSNP, c("BFMI861.S1", "BFMI861.S1.1")], chr[dSNP, nonDiabetic])

    write.table(snpOUT, file=paste0("MM_diffences_chr", x, "_snps.txt"), sep="\t",row.names = FALSE)            # Write individual output files
    
    if(x==1){
      write.table(snpOUT, file=combinedFN, sep="\t",row.names = FALSE)                                          # Write the combined output files *with header*
    }else{
      write.table(snpOUT, file=combinedFN, sep="\t",row.names = FALSE, append=TRUE, col.names=FALSE)            # Write the combined output files *no header, and appending*
    }
    cat(paste0("Done with chromosome ", x, ", detected: ", length(dSNP), " differences\n"))
  }else{
    cat("Skipping chromosome", x, "no consistent SNPs\n")
  }
})

# Might be interesting as targets, build 38 coordinates :X
#
# Insr: Insuline receptor       Chromosome 8:    3,150,922 -   3,279,617 reverse strand.
# GLUT4/Slc2a4 :                Chromosome 11:  69,942,539 -  69,948,188
#
# From the article: Stratigopoulos et al. (Cell metabolism)
# 
# FTO:                          Chromosome 8:   91,313,532 -  91,668,439
# Rpgrip1l:                     Chromosome 8:   91,217,030 -  91,313,262 reverse strand
# Cux1:                         Chromosome 5:  136,248,135 - 136,567,490 reverse strand
# Lepr:                         Chromosome 4:  101,717,404 - 101,815,352
# STAT3:                        Chromosome 11: 100,885,098 - 100,939,540 reverse strand
#
# Npy:                          Chromosome 6:   49,822,710 -  49,829,507
# Agrp:                         Chromosome 8:  105,566,700 - 105,568,298 reverse strand
# Pomc:                         Chromosome 12:   3,954,951 -   3,960,618
# Mc4r:                         Chromosome 18:  66,857,715 -  66,860,472 reverse strand

interestingTargets <- rbind(c("Insr",          8,   3150922,   3279617),
                            c("GLUT4/Slc2a4", 11,  69942539,  69948188),
                            c("FTO",           8,  91313532,  91668439),
                            c("Rpgrip1l",      8,  91217030,  91313262),
                            c("Cux1",          5, 136248135, 136567490),
                            c("Lepr",          4, 101717404, 101815352),
                            c("STAT3",        11, 100885098, 100939540),
                            c("Npy",           6,  49822710,  49829507),
                            c("Agrp",          8, 105566700, 105568298),
                            c("Pomc",         12,   3954951,   3960618),
                            c("Mc4r",         18,  66857715,  66860472))

S1versusBFMI <- read.table("S-1versusBFMI.txt", header=TRUE)

for(tid in 1:nrow(interestingTargets)){
  target   <- interestingTargets[tid,]
  chrSNPs  <- which(S1versusBFMI[,"Chromosome"] == as.numeric(target[2]))
  snpOnCHR <- S1versusBFMI[chrSNPs, ]
  inRange  <- ( snpOnCHR[, "Pos.mm10"] > as.numeric(target[3])- 20000 & snpOnCHR[, "Pos.mm10"] < as.numeric(target[4]) + 200000 )
  if(any(inRange)){
    snpsInGene <- snpOnCHR[which(inRange), ]
    cat("Found a SNP in:",target[1],", SNP:", as.character(snpsInGene[,"dbSNP.RS.ID"]),"\n")
    return()
  }
}

# Inside the gene regions:
#
# Found a SNP in: Cux1 , SNP: rs32248143  5:   136051558    -> Upstream gene variant located in 5'
# Found a SNP in: Pomc , SNP: rs49073952 12:     3957608    -> intron variant in Pomc & Downstream gene variant in EFR3

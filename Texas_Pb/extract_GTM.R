setwd("D:/Edrive/Mouse/Texas_Pb")

gt1 <- read.csv("input/Texas_A&M_Threadgill_MURCOMV02_20191118_FinalReport.txt", header=TRUE, skip = 9, sep="\t", na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
gt2 <- read.csv("input/Texas_A&M_Threadgill_MURCOMV02_20200519_FinalReport.txt", header=TRUE, skip = 9, sep="\t", na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
gt3 <- read.csv("input/Texas_A&M_Threadgill_MURCOMV02_20200526_FinalReport.txt", header=TRUE, skip = 9, sep="\t", na.strings=c("","-", "na", "NA", "NaN", "X", "x"))

gts <- rbind(gt1,gt2,gt3)

markers <- unique(gts[,1])
samples <- unique(gts[,2])

# Filter SNPs based on 0.7, as recommended by illumina:
# - https://www.illumina.com/Documents/products/technotes/technote_gencall_data_analysis_software.pdf

gtm <- matrix(NA, length(markers), length(samples), dimnames=list(markers, samples))

for(x in 1:nrow(gts)){
  if(!is.na(gts[x, "GC.Score"]) && gts[x, "GC.Score"] > 0.7){
    gtm[gts[x, 1],gts[x, 2]] <- paste0(gts[x, "Allele1...Forward"], gts[x, "Allele2...Forward"])
  }
}

write.table(gtm, "genotypes_all.txt", sep = "\t", quote = FALSE)

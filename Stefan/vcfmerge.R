# vcfmerge.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Functions to merge novelSNPer output, with known RS ids from vcf files

# Try to correct VCF towards the novelSNPer reference
matchref <- function(current.data, refallele, rsnumber){
  if(current.data[4] == refallele) return(paste0(current.data[4], "/", current.data[5]))
  if(current.data[5] == refallele) return(paste0(current.data[5], "/", current.data[4]))
  cat("WARN: Ref allele does not match for", rsnumber, ", ref:", refallele, ", vcf data:",current.data[4], "and", current.data[5], "\n")
  return(paste0(current.data[4],"/",current.data[5]))
}

# Get the real RS positions, and not use the weird ones from the VCF file
getRSPOS <- function(column8){ return(as.integer(unlist(strsplit( (unlist(strsplit(column8, split = ";")))[1], split = "="))[2])) }

chrs <- c(seq(19,2,-1), "X")

for(chr in chrs){
  novelsnper  <- read.csv(paste0("/home/neubert/Mouse_NGS_MDC_860v2/mm10/NovelSNPer/860v2_chr", chr, "_NovelSNPer_detailed.out"), sep="\t")   # Output form novel SNPer
  vcffile     <- paste0("/home/neubert/ref_genomes/mm10/Mus_musculus_GRCm38_dbSNP137_snpome/vcfchr", chr, ".vcf")

  novelout <- cbind("", TRUE, "", novelsnper)                                             # The output matrix
  colnames(novelout) <- c("RS", "isNovel", "VCF", colnames(novelsnper))                   # We're adding 3 columns to it
  novelout <- as.matrix(novelout)                                                         # everything as.character, so we don't get factors involved

  con <- file(vcffile) 
  open(con)
  current.line <- 1

  while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    current.data  <- unlist(strsplit(line, split = "\t"))
    location      <- getRSPOS(current.data[8])
    if(location %in% novelsnper[, "Start"]){
      locNovelSNP <- which(novelsnper[, "Start"] == location)

      refallele <- as.character(novelsnper[locNovelSNP[1], "ref_Allele"])
      # cat("\nMatching location found:", location, "at line:", current.line, "at novelSNP:", locNovelSNP, ", RS:", current.data[3], "\n")

      vcfallele <- matchref(current.data, refallele, current.data[3])
      vcfallele <- gsub(",", "/", vcfallele)
      # cat(" - VCF allele:", vcfallele, "\n")
      for(lNovelSNP in locNovelSNP){
        novelSNP <- paste0(refallele, "/", as.character(novelsnper[lNovelSNP, "Allele1"]))
        # cat(" - NovelSNP allele:", novelSNP, "\n")
        novelout[lNovelSNP, 1] = as.character(current.data[3])
        novelout[lNovelSNP, 2] = (novelSNP != vcfallele)
        novelout[lNovelSNP, 3] = vcfallele
      }
    }
    current.line  <- current.line + 1
  }
  close(con)

  write.table(novelout, paste0("860v2_chr", chr, "_NovelSNPer_detailed_RS.out"), sep="\t", row.names=FALSE)   # Output form novel SNPer
}

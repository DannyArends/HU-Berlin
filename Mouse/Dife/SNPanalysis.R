# Analysis of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2016
# first written Jun, 2016

setwd("~/NAS/Mouse/DNA/Sequencing/DifeMouse")

samples <- c("ext_L7254", "ext_L7256", "ext_L7257", "ext_L7255", "ext_L7258")
names(samples) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")

line.n  <- 1
Tfile   <- file("ALL_variants.raw.joint.g.vcf.snpeff_annotated_GRCm38.79_FILTERED.vcf", "r")                    # Open a connection
hasHeader <- FALSE
while(length((line = readLines(Tfile, n = 1))) > 0) {
   # cat(line, "\n")
    
    if(!hasHeader && substr(line, 1, 2) == "#C"){
      #cat("header at", line.n, "\n")
      header <- strsplit(substr(line, 2, nchar(line)), "\t")[[1]]
      hasHeader <- TRUE
      headerInSample <- which(header %in% samples)
      cat(paste0(c("ID", "Chr", "Pos", "ref", "alt", "q", header[which(header %in% samples)]), sep="", collapse="\t"), "\n", sep="", file="out.txt")
    }else if(hasHeader){
      #cat("We have the header at", line.n, "\n")
      allrowdata <- strsplit(line, "\t")[[1]]
      cleanedGeno <- substr(allrowdata[headerInSample], 1, 3)
      cat(paste0(c(allrowdata[c(3,1,2,4,5,6)], cleanedGeno), sep="", collapse="\t"), "\n",sep="", file="out.txt", append = TRUE)
    }
    line.n <- line.n + 1
  if(line.n %% 10000 == 0) cat("at ", line.n, "\n")
}
close(Tfile)


# VCFcheck.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Small script to test SNP location data from a .vcf file

setwd("D:/Stefan_Mouse_F3/VCF")

chr <- 1
vcffile <- paste0("vcfchr", chr, ".vcf")

getRSPOS <- function(column8){ return(as.integer(unlist(strsplit( (unlist(strsplit(column8,";")))[1], "="))[2])) }

con <- file(vcffile) 
open(con);
current.line <- 1
while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  current.data  <- unlist(strsplit(line, split="\t"))
  location      <- as.integer(current.data[2])
  loc2          <- getRSPOS(current.data[8])
  if(location != loc2){
    cat(current.line,"",RSPOSidx, ": location in VCF is Borked", location,"!=", loc2,"\n")
    #return()
  }
  current.line  <- current.line + 1
}
close(con)

### Selecting SNPs in our chromosome 3 regions of interest 


#setwd("E:/Mouse/DNA/Keane/RawData")

cat("", file="/home/arends/output.txt")

line.n  <- 1
Tfile <- file("v5redSNPheaderChr3.vcf", "r")
while(length((line = readLines(Tfile, n = 1) )) > 0){                         # Read a line, if available
  if((line.n %% 10000) == 0)cat("line:", line.n, "\n")
  if(substr(line,1,1) != "#"){
    elements <- strsplit(line, "\t")[[1]]
    loc <- as.numeric(elements[2])
    if(is.na(loc)){
      geno <- unlist(lapply(strsplit(elements[10:length(elements)], ":"), "[", 1))
      cat(c(elements[1:5], geno), "\n", file="/home/arends/output.txt",append=TRUE)
    }else if((loc > 14841429 && loc < 17563072) || (loc > 35308230 && loc < 36854743)){
      genoC <- unlist(strsplit(paste0(elements[c(4,5)],collapse=","), ","))
      geno <- unlist(lapply(strsplit(elements[10:length(elements)], ":"), "[",1))

      if(all(grepl("/", geno))){
        genoChar <- unlist(lapply(lapply(lapply(strsplit(geno, "/"),as.numeric), "+", 1), function(x){
          if(length(x) > 1){
            if(!is.na(x[1]) && !is.na(x[2])){
              if(x[1] != x[2]) return("H")
            }
          }
          return(genoC[x][1])
        }))
        cat(c(elements[1:5], genoChar), "\n", file="/home/arends/output.txt",append=TRUE)
      }else{
        cat("DEFNKT, SHOULD NOT GET HERE\n")
      }
    }
  }
  line.n <- line.n + 1
}

genotypes <- read.table("genotypesRegionsChr3.txt", header=TRUE, sep=" ")
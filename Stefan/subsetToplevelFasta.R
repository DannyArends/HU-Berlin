### Selecting SNPs in our chromosome 3 regions of interest 

chromosomes <- c(1:19,"X","Y","MT")

cat("", file="/home/arends/RNASeq/Reference/Mus_musculus.GRCm38.dna.fa")
copy <- FALSE
line.n  <- 1
Tfile <- gzfile("/home/arends/RNASeq/Reference/Mus_musculus.GRCm38.dna.toplevel.fa.gz", "r")
while(length((line = readLines(Tfile, n = 1) )) > 0){                         # Read a line, if available
  #if((line.n %% 100000) == 0) cat("line:", line.n, "\n")
  if(substr(line,1,1) == ">"){
    chr <-  gsub(">", "", strsplit(line," ")[[1]][1])
    if(chr %in% chromosomes){
      copy <- TRUE
    }else{
      copy <- FALSE
    }
    cat(line, copy,"\n")
  }
  if(copy) cat(paste0(line, "\n", collapse=""), sep="", file="/home/arends/RNASeq/Reference/Mus_musculus.GRCm38.dna.fa", append=TRUE)
  line.n <- line.n + 1
}

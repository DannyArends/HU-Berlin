# miRNA analysis on Cows 3' UTR
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2015
# first written Aug, 2015

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

execute("wget http://genie.weizmann.ac.il/pubs/mir07/64bit_exe_pita_prediction.tar.gz", "64bit_exe_pita_prediction.tar.gz")
execute("tar xvfz 64bit_exe_pita_prediction.tar.gz", "Makefile")
execute("make install", "pita_prediction.pl")

execute(paste0("./pita_prediction.pl -utr ", file.UTR, "-mir ", file.miRNA," -prefix ", outputname))
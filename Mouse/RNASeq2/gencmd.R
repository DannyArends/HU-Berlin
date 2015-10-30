# Generate commands for the RNAseq pipeline, for all the samples in sampledescription.txt
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sept, 2015
# first written Juli, 2015

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/Experiment")
samples <- read.table("sampledescription.txt", sep="\t", header=TRUE, colClasses="character")

fdone <- c("4422","4423","4424","4425","4426","4427","4868","5067","5068","5069","5070","5074","5075","5076","5077","5078","5079","5080","5081","6339","6340","6341","6342","6343","6344","6345","6346")
done  <- c("4422","4426","5080","4424","5081","5079","5076","5078","5077","5075","5072","5070","5068","5071","6345","6339","6341","6343","6344","5073","4427","6346","6342","5074")

for(x in 1:nrow(samples)){
  fname <- paste0(samples[x,"Lib_id"], "_", samples[x,"TagLane"])
  if(!samples[x,"Lib_id"] %in% fdone){
    cat(paste0("nohup Rscript pipeline.R /home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/FASTQ/", fname," > ",fname,".log 2>&1&"),"\n")
  }
}

done[!done %in% fdone]  # Have stopped somewhere 1/2 way


"5073" "5075" "5072" "5070" "5071" "5068" "4425" "4423" # Have stopped somewhere 1/2 way
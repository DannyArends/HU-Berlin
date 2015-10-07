# Generate commands for the RNAseq pipeline, for all the samples in sampledescription.txt
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sept, 2015
# first written Juli, 2015

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/Experiment")
samples <- read.table("sampledescription.txt", sep="\t", header=TRUE, colClasses="character")
done <- c("4422","4423"  ,"5080"  ,"4424"  ,"5081"  ,"5079"  ,"5076"  ,"5078"  ,"5077"  ,"5075"  ,"5072"  ,"5070"  ,"5068"  ,"5071"  ,"6345"  ,"6339"  ,"6341"  ,"6343"  ,"6344"  ,"5073"  ,"4427"  ,"6346"  ,"6342" )

for(x in 1:nrow(samples)){
  fname <- paste0(samples[x,"Lib_id"], "_", samples[x,"TagLane"])
  if(!samples[x,"Lib_id"] %in% done){
    cat(paste0("nohup Rscript pipeline.R /home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/FASTQ/", fname," > ",fname,".log 2>&1&"),"\n")
  }
}

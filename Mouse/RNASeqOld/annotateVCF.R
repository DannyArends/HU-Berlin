#
# Annotate a VCF using GTF information
#


setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
line.n <- 1
file.GTF          <- file("GTF/Mus_musculus.GRCm38.81.gtf", "r")               # Open the GTF file connection
file.cannonical   <- read.table("GTF/mouse_canonical_transcripts_v81.fa", sep="\t", check.names=FALSE, colClasses="character")

if(!file.exists("GTF/ExonMatrix.GRCm38.81.txt")){
  file.output       <- file("GTF/ExonMatrix.GRCm38.81.txt", "w")                 # Open the output connection
  cnt.nC <- 0 ; cnt.nE <- 0                                                      # Some statistics

  cat(paste("chr","start","stop","gene_id","exon_id","transcript_id","exon_number","gene_name",sep="\t"),"\n",sep = "", file = file.output)
  while(length((line = readLines(file.GTF, n = 1))) > 0){ # Read a single line cat(line, "\n")
    if(substr(line,1,1) != "#"){
      line.split <- strsplit(line, "\t")[[1]]
      if(line.split[3] == "exon"){
        transcript_id <- sub(".*?transcript_id \"(.*?)\";.*", "\\1", line.split[9])
        if(transcript_id %in% file.cannonical[, 2]){
          cat(paste(line.split[1], line.split[4], line.split[5],
                    sub(".*?gene_id \"(.*?)\";.*", "\\1", line.split[9]), 
                    sub(".*?exon_id \"(.*?)\";.*", "\\1", line.split[9]), 
                    transcript_id, 
                    sub(".*?exon_number \"(.*?)\";.*", "\\1", line.split[9]), 
                    sub(".*?gene_name \"(.*?)\";.*", "\\1", line.split[9]), sep = "\t"),"\n", sep="", file=file.output, append=TRUE)
        }else{
        cnt.nC <- cnt.nC + 1
        }
      }else{
        cnt.nE <- cnt.nE + 1
      }
    }
    line.n <- line.n + 1
    if(line.n %% 1000 == 0) cat("Done", line.n, "lines",cnt.nE,"non-exons", cnt.nC,"non-canonical\n")
  }
}else{
  cat("Skipping building the file\n");
}
close(file.GTF)
close(file.output)
exon.matrix <- read.table("GTF/ExonMatrix.GRCm38.81.txt", sep="\t", header = TRUE)
exon.matrix[1:5,]

file.VCF <- read.table("ReAnalysisSNPs/population.vcf", header = TRUE, colClasses="character")

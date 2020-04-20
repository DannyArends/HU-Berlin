cmdlineargs <- commandArgs(trailingOnly = TRUE)
seqID <- as.character(cmdlineargs[1])          # sequenceID
sampleID <- as.character(cmdlineargs[2])       # sampleID
R1 <- as.character(cmdlineargs[3])             # fastq R1 file
R2 <- as.character(cmdlineargs[4])             # fastq R2 file

outputfolder <- paste0("/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/", seqID , "/")
referencefolder <- "/home/arends/NAS/Mouse/Reference_Genomes/GRCm38_99/"
reference <- paste0(referencefolder, "/Mus_musculus.GRCm38.dna.toplevel.fa.gz")
logfile <- paste0(outputfolder, "log.txt")

cat("--------------------------------------------------------------------\n")
cat("seqID:", seqID, "\n")
cat("sampleID:", sampleID, "\n")
cat("R1:", R1, "\n")
cat("R2:", R2, "\n")
cat("outputfolder:", outputfolder, "\n")
cat("logfile:", logfile, "\n")
cat("reference:", reference, "\n")
cat("--------------------------------------------------------------------\n")

dir.create(outputfolder, showWarnings = FALSE, recursive = TRUE)
cat("", file = logfile)
cat("Created output folder\n", file = logfile, append = TRUE)

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

### Trimmomatic: Remove adapters and trim reads based on quality scores (1 to 2 hours) ###
trimmomatic <- "/home/arends/Github/trimmomatic/trimmomatic-0.38/dist/jar/"
adapters <- "/home/arends/Github/trimmomatic/trimmomatic-0.38/adapters/TruSeq3-PE.fa"

trimmomaticout <- paste0(outputfolder, "trimmed/")
dir.create(trimmomaticout, showWarnings = FALSE, recursive = TRUE)
cat("Created Trimmomatic output folder\n", file = logfile, append = TRUE)

outputfiles <- c(paste0(trimmomaticout, sampleID, "1.P_trimmed.fastq.gz"), paste0(trimmomaticout, sampleID, "1.U_trimmed.fastq.gz"), 
                 paste0(trimmomaticout, sampleID, "2.P_trimmed.fastq.gz"), paste0(trimmomaticout, sampleID, "2.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.38.jar PE")
params  <- paste0("ILLUMINACLIP:", adapters, ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

command <- paste(cmdBase, R1, R2, outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)

cat(command, "\n", file = logfile, append = TRUE)

q("n")



setwd("~/NAS/Goat/DNA/Sequencing/fastq")
snames <- read.table("../sample_names.txt", header=TRUE, colClasses="character")

reads1 <- gzfile("lane7_R1.fastq.gz", "r")
barcode1 <- file("lane7_R2.fastq.gz", "r")
barcode2 <- gzfile("lane7_R3.fastq.gz", "r")
reads2 <- gzfile("lane7_R4.fastq.gz", "r")

cnt <- 1

for(fstem in snames[, 1]){
  cat("", file = paste0("DeMux/", fstem, "_R1.fastq"))
  cat("", file = paste0("DeMux/", fstem, "_R2.fastq"))
}

ours <- 0
unkn <- 0
barcodes <- c()

while ( TRUE ) {
  bc1 <- readLines(barcode1, n = 4)
  bc2 <- readLines(barcode2, n = 4)
  read1 <- readLines(reads1, n = 4)
  read2 <- readLines(reads2, n = 4)
  if ( length(bc1) == 0 ) {
    break
  }
  inS <- which(snames[,2] == bc1[2] & snames[,3] == bc2[2])
  if(cnt %% 10000 == 0) cat("Done", paste0(bc1[2],"-",bc2[2]), ours, "+", unkn, " = ", cnt, "Reads = ", inS, "\n")
  if(length(inS) == 1) {
    fstem <- paste0("DeMux/", snames[inS, 1])
    cat(paste0(paste(read1, collapse="\n"),"\n"), file = paste0(fstem, "_R1.fastq"), append=TRUE)
    cat(paste0(paste(read2, collapse="\n"),"\n"), file = paste0(fstem, "_R2.fastq"), append=TRUE)
    ours <- ours + 1
  } else {
    unkn <- unkn + 1
  }
  cnt <- cnt + 1
}
q("no")

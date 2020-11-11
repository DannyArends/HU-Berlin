setwd("")
mHol <- read.csv("chr25_Holstein.frq", sep="\t", header=FALSE, skip=1, row.names=NULL)
mHolF <- read.csv("chr25_HolsteinFriesian.frq", sep="\t", header=FALSE, skip=1, row.names=NULL)
mGelb <- read.csv("chr25_Gelbvieh.frq", sep="\t", header=FALSE, skip=1, row.names=NULL)
mBraun <- read.csv("chr25_OriginalBraunvieh.frq", sep="\t", header=FALSE, skip=1, row.names=NULL)

allB <- rbind(mHol, mHolF, mGelb, mBraun)
positions <- unique(allB[,2])

#2b) Reading a text file line by line
line.n  <- 1
nunique <- 0
chrPos <- c()
Tfile <- gzfile("/home/danny/NAS/Cattle/DNA/Variants/DSN/raw/Chr25/raw.vcf.gz", "r") # Get a pointer to the file
while(length((line = readLines(Tfile, n = 1) )) > 0){ # Read a line, if available
  splitted <- strsplit(line,"\t")[[1]]
  line.length <- length(splitted)
  if (line.length == 314) {
    #cat("line:", line.n, "length:", line.length, "\n")                                # Number of words per line
    pos <- splitted[2]
    if(!pos %in% positions){
      nunique <- nunique + 1
      chrPos <- c(chrPos, pos)
    }
  }
  line.n <- line.n + 1                                                                                  # Increase out line number counter
}
close(Tfile)


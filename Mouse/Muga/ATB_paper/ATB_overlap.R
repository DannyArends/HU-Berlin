
setwd("D:/Edrive/Mouse/DNA/MegaMuga/")

patRegions <- read.table("Analysis/ATB_PAT.txt",sep="\t",header=TRUE, check.names=FALSE)
matRegions <- read.table("Analysis/ATB_MAT.txt",sep="\t",header=TRUE, check.names=FALSE)

patUnique <- 0
patOverlaps <- c()
for(x in 1:nrow(patRegions)){
  chr <- patRegions[x, "Chr"]
  spos <- patRegions[x, "Start"]
  epos <- patRegions[x, "Stop"]
  onChr <- which(matRegions[,"Chr"] == chr)
  overlap <- FALSE
  owith <- NULL
  for (y in onChr) {
    if(matRegions[y, "Start"] >= spos && matRegions[y, "Start"] <= epos){
      overlap <- TRUE
      owith <- y
    }
    if(matRegions[y, "Stop"] >= spos && matRegions[y, "Stop"] <= epos){
      overlap <- TRUE
      owith <- y
    }
  }
  if(!overlap){
     patUnique <- patUnique + 1
  }else{
    patOverlaps <- rbind(patOverlaps, c(x, owith))
  }
}

matUnique <- 0
matOverlaps <- c()
for(x in 1:nrow(matRegions)){
  chr <- matRegions[x, "Chr"]
  spos <- matRegions[x, "Start"]
  epos <- matRegions[x, "Stop"]
  onChr <- which(patRegions[,"Chr"] == chr)
  overlap <- FALSE
  owith <- NULL
  for (y in onChr) {
    if(patRegions[y, "Start"] >= spos && patRegions[y, "Start"] <= epos){
      overlap <- TRUE
      owith <- y
    }
    if(patRegions[y, "Stop"] >= spos && patRegions[y, "Stop"] <= epos){
      overlap <- TRUE
      owith <- y
    }
  }
  if(!overlap){
     matUnique <- matUnique + 1
  }else{
    matOverlaps <- rbind(matOverlaps, c(x, owith))
  }
}

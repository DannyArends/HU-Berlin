#
# New QTL mapping of BFMI, using a slightly different population structure correction
#

setwd("E:/Mouse/DNA/MegaMuga")

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genoF2 <- read.table("extra/genotypesF2.txt", sep="\t", check.names=FALSE, colClasses="character")                  # Normal A, H, B genotypes
genoF1 <- read.table("extra/genotypesF1.txt", sep="\t", check.names=FALSE, colClasses="character")                  # Normal A, H, B genotypes
phenotypes <- read.table("extra/phenotypes.txt", sep="\t", check.names=FALSE, colClasses="character")               # Phenotypes

F2  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                       # The F2 individuals
F1  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]                                                      # The F1 individuals

d63 <- as.numeric(phenotypes[F2,"d63"])

todo <- which(!(rownames(genoF2) %in% done))
todo <- todo[-length(todo)]
from = 20000; to = 21259; fname = "extra/ALL_MxM_interactionsTODO3.txt"

## Scan interactions between m1 and the end, taken from a todo list
cat(paste0(rownames(genoF2),collapse="\t"), "\n", file=fname, sep="")
for(m1i in todo) {
  m2i <- m1i + 1
  m1data <- as.factor(unlist(genoF2[m1i, F2]))
  res <- rep(NA, 21260)
  for(m2i in m1i:21260) {
    tryCatch( {
        model <- lm(d63 ~ m1data * as.factor(unlist(genoF2[m2i, F2])))
        aovm <- anova(model)[[5]]
        if(length(aovm) == 4) res[m2i] <- aovm[3]
      }, error = function(e) {
        cat(paste("Error caught:", e, "\n"))
      }
    )
  }
  cat(paste0(rownames(genoF2)[m1i], "\t", paste0(res, collapse="\t"), collapse=""), "\n", file=fname, sep="", append=TRUE)
  cat(m1i, "Done\n")
  m1i <<- m1i + 1
}

## See what finished, by scannig the rownames of the files... used to split up work
setwd("E:/Mouse/DNA/MegaMuga")
intfiles <- dir("extra")[grep("ALL", dir("extra"))]
done <- c()
for(mf in intfiles){
  cat(mf,"\n")
  #2) Reading a text file line by line
  line.n  <- 1
  maxlod <- 0
  Tfile <- file(paste0("extra/",mf), "r")
  while(length((line = readLines(Tfile, n = 1) )) > 0){                                                   # Read a line, if available
    if(line.n > 1){
      splitted <- strsplit(line,"\t")[[1]]
      #cat("line:", line.n, "length:", length(splitted), "\n")                                # Number of words per line
      values <- -log10(as.numeric(splitted[-1]))
      mlod <- max(values,na.rm=TRUE)
      if(mlod > maxlod){
        #break;
        idx <- which(values == mlod)
        cat(splitted[1], "(", which(header == splitted[1]), ") ", header[idx], " - ", mlod, "\n")
        maxlod <- mlod
      }
      done <- c(done, splitted[1])
    }else{
      header <- strsplit(line,"\t")[[1]]
    }
    line.n <- line.n + 1
  }
}

## Load the created pvalue files, transform to LOD and filter
setwd("E:/Mouse/DNA/MegaMuga")
intfiles <- dir("extra")[grep("ALL", dir("extra"))]
outfile <- "extra/LODinteractions.txt"
cat("", file = outfile)
noHeader <- TRUE
for(mf in intfiles){
  cat(mf,"\n")
  #2) Reading a text file line by line
  line.n  <- 1
  Tfile <- file(paste0("extra/",mf), "r")
  while(length((line = readLines(Tfile, n = 1) )) > 0){                                                   # Read a line, if available
    if(line.n > 1){
      splitted <- strsplit(line,"\t")[[1]]
      #cat("line:", line.n, "length:", length(splitted), "\n")                                # Number of words per line
      values <- -log10(as.numeric(splitted[-1]))
      mlod <- max(values,na.rm=TRUE)
      if(mlod >= 5) {
        cat(splitted[1],"\t", paste0(round(values,1), collapse="\t"), "\n",sep="", file = outfile, append=TRUE)
      }
    }else{
      if(noHeader) {
        cat(line, "\n", file = outfile, append=TRUE)
        header <- 
        noHeader <- FALSE
      }
    }
    line.n <- line.n + 1
  }
}

## Load in the higest lod scores, and do some more filtering
setwd("E:/Mouse/DNA/MegaMuga")
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
highlods <- read.csv("extra/LODinteractions.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
nNA <- apply(highlods, 2, function(x){sum(is.na(x))})
highlods <- highlods[,-which(nNA == nrow(highlods))]

above7r <- apply(highlods, 1, function(x){ any(x >= 7,na.rm=TRUE); })
above7c <- apply(highlods, 2, function(x){ any(x >= 7,na.rm=TRUE); })

signsubset <- apply(highlods[above7r, above7c], 2, as.numeric)
rownames(signsubset) <- rownames(highlods[above7r, above7c])

rO <- sort(unlist(lapply(rownames(signsubset), function(x){which(rownames(map) == x)})),i=TRUE)$ix
cO <- sort(unlist(lapply(colnames(signsubset), function(x){which(rownames(map) == x)})),i=TRUE)$ix

chrsR <- (which(diff(as.numeric(map[rownames(signsubset[rO,cO]),c(1)])) > 0)-0.5) / (length(rO)-1)
chrsC <- (which(diff(as.numeric(map[colnames(signsubset[rO,cO]),c(1)])) > 0)-0.5) / (length(cO)-1)

image(signsubset[rO,cO], breaks = c(0, 3, 6, 9, 12, 15), col=c("white", "gray", "yellow", "orange", "red"), xaxt='n', yaxt='n')
axis(1, seq(0,1,1 / (length(rO)-1)), paste0(map[rownames(signsubset[rO,cO]),c(1)], "-", round(as.numeric(map[rownames(signsubset[rO,cO]),c(2)]) / 1000000,1)), las=2,cex.axis = 0.9)
axis(2, seq(0,1,1 / (length(cO)-1)), paste0(map[colnames(signsubset[rO,cO]),c(1)], "-", round(as.numeric(map[colnames(signsubset[rO,cO]),c(2)]) / 1000000,1)), las=2,cex.axis = 0.7)
abline(v=chrsR)
abline(h=chrsC)
box()


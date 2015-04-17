# Analysis of the microarray data from Agilent
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

arrays <- c("File1.txt", "File2.txt")     # file names

rawdata <- NULL
for(filename in arrays){
  marray <- read.csv(paste0("Location/", filename), sep = "\t", skip = 9, header = TRUE)      # I assume files are in workdir/Location/File1.txt
  rawdata <- cbind(rawdata, marray$gProcessedSignal)                                          # Get the processed intensity signals
}
rawdata <- cbind(marray[,c("ProbeName", "Sequence")], rawdata)
colnames(rawdata) <- c("ProbeName", "Sequence", arrays)                                       # Add the name and the sequence of the probe

#Some probes are on the array multiple times, so let's de-duplicate them (we take the average expression over the similar probes)
uniqueprobes <- unique(rawdata[,"ProbeName"])
arraydata <- matrix(NA, length(uniqueprobes), (ncol(rawdata)-1), dimnames = list(uniqueprobes, colnames(rawdata)[-1]))
cat("De-duplicating probes\n")
for(probe in uniqueprobes){
  prows <- which(rawdata[,"ProbeName"] == probe)
  probeseq <- rawdata[prows,"Sequence"][1]
  arraydata[probe, ] <- c(as.character(probeseq), apply(rawdata[which(rawdata[,"ProbeName"] == probe), arrays], 2, mean))
}
arraydata <- data.frame(arraydata)
write.table(arraydata, file="Analysis/raw_arraydata.txt", sep="\t", row.names=TRUE)

# Preprocessing
arraydata[,arrays] <- log2(arraydata[,arrays])                              # Log2 transformation
boxplot(arraydata[,arrays], las=2)
arraydata[,arrays] <- normalize.quantiles(as.matrix(arraydata[,arrays]))    # Quantile normalisation
boxplot(arraydata[,arrays], las=2)
write.table(arraydata, file="Analysis/normalized_log2_arraydata.txt", sep="\t", row.names=TRUE)

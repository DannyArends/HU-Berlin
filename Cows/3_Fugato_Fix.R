
## Fix the wrong annotation for fugato (120 animals are missing in action)
sampleInfo <- read.csv("Fugato/sampleinfoKielAna.txt",sep="\t", header=TRUE, colClasses="character", stringsAsFactor=FALSE, check.names=FALSE)
fugato <- read.table("Fugato/genotypes.txt",check.names=FALSE, colClasses='character',stringsAsFactor=FALSE)
write.table(sampleInfo[colnames(fugato),], file="Fugato/sampleInfo_matched.txt", sep="\t", quote=FALSE)


# Analysis of JAX known/pseudo gene data to create a custom Agilent array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#

setwd("E:/Mouse/RNA/ArrayDesign/Agilent")

genedata <- read.table("MGI_Gene_Model_Coord.rpt", sep="\t", header=TRUE, row.names=1)                    # All genes in the mouse genome according to MGI
genesOnChr3 <- which(genedata[,"Ensembl.gene.chromosome"]==3)                                             # Take the ones on chromosome 3
genedataChr3 <- genedata[genesOnChr3,]

write.table(rownames(genedataChr3),"genesOnChr3.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

cat("Number of known/predicted genes:", length(genesOnChr3),"\n")

# We can now upload this to MGI batch query: http://www.informatics.jax.org/batch
# For Agilent, we only need the: GenBank/RefSeq ID of the gene

genbankFromMGI <- read.table("MGIBatchReport_20140701_032339.txt", sep="\t", header=TRUE)

knownrna <- grep("NM_", genbankFromMGI[,4])
predirna <- grep("XM_", genbankFromMGI[,4])

cat("Known mRNA transcripts :", length(knownrna),"\n")
cat("Predicted mRNA transcripts :", length(predirna),"\n")

rnaToTarget <- c(as.character(genbankFromMGI[knownrna,"GenBank.RefSeq.ID"]), as.character(genbankFromMGI[predirna,"GenBank.RefSeq.ID"]))
write.table(rnaToTarget ,"rnasForAgilenteArray.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

# We can now upload this to Agilent eArray: https://earray.chem.agilent.com/
# Select GE Probe Design, -> Genbank Accession & M. Musculus (UniGene Build #186)

# Which ones can we leave off the array
setwd("E:/Mouse/RNA/ArrayDesign/Agilent")

onArray <- read.csv("PGRID172954COMPLETE.tdt", sep="\t", header=TRUE, colClasses=c("character"))
notExpressed <- read.table("MGItoREFSEQ.txt", sep="\t", header=TRUE, colClasses=c("character"))

targetOnArray <- na.omit(unlist(lapply(strsplit(onArray[,"TargetID"],"|",fixed=TRUE),"[",2)))
removalCandidates <- targetOnArray[targetOnArray %in% notExpressed[,"GenBank.RefSeq.ID"]]

toRem <- onArray[which(onArray[,"TargetID"] %in% paste0("ref|", removalCandidates)), ]
toRem <- toRem[-grep("chr3", toRem[,"ChromosomalLocation"]),]

# Which ones are already on the array

addToArray <- read.csv("OurProbeGroup.tdt", sep="\t", header=TRUE, colClasses=c("character"))
wantToAdd <- unlist(lapply(strsplit(na.omit(unlist(lapply(strsplit(addToArray[,"TargetID"],"|",fixed=TRUE),"[",2))),".", fixed=TRUE),"[",1))

toAdd <- wantToAdd[which(!wantToAdd %in% targetOnArray)]
cat("We want to add", length(toAdd), "probes / genes on chromosome 3\n")

# We find: 1457 probes not on the array yet, which we want to add
# We find: 2074 probes that can be removed/swapped by our new probes

# Create an output TDT file containing our probes, and the ones from Agilent
set.seed(1)
doNotTake <- which(onArray[,"TargetID"] %in% toRem[,"TargetID"])
cat("We could remove", length(doNotTake), "non expressed genes\n")
notTaken <- sample(doNotTake, 1457)
AgilentPartRemoved <- onArray[notTaken,]                                                  # Because of duplicates we only need 1457 ID's
AgilentPart <- onArray[-notTaken,]                                                        # Because of duplicates we only need 1457 ID's
HUBerlinPart <- addToArray[which(wantToAdd %in% toAdd),]

refseqIDs <- unlist(lapply(strsplit(na.omit(unlist(lapply(strsplit(HUBerlinPart[,"TargetID"],"|",fixed=TRUE),"[",2))),".", fixed=TRUE),"[",1))
MGIs <- genbankFromMGI[match(refseqIDs,genbankFromMGI[,"GenBank.RefSeq.ID"]),"MGI.Gene.Marker.ID"]

HUBerlinPart[,"GeneSymbol"] <- genedataChr3[MGIs,"marker.symbol"]
HUBerlinPart[,"Description"] <- genedataChr3[MGIs,"marker.name"]
HUBerlinPart[,"ChromosomalLocation"] <- paste0("mm|chr", genedataChr3[MGIs,"Ensembl.gene.chromosome"], ":", genedataChr3[MGIs,"Ensembl.gene.start"], "-", genedataChr3[MGIs,"Ensembl.gene.end"])

newProbeGroup <- rbind(AgilentPart, HUBerlinPart)

IDS <- c(2086, 2385, 3241, 5977, 18361, 18973, 19595, 25544, 30197, 33810, 37288, 37944, 43889, 44526, 47079, 48497, 50800, 51305)-1
newProbeGroup[IDS, "Description"] <- "Agilent eArray bug for this description"

write.table(newProbeGroup, "HUBerlinProbeGroup.tdt", sep="\t", quote=FALSE, row.names=FALSE)
dim(AgilentPart)
dim(HUBerlinPart)
dim(onArray)
dim(newProbeGroup)


write.table(onArray[doNotTake, ],"PossibleCandidatesForRemoval.txt",sep="\t", row.names=FALSE)
write.table(onArray[notTaken, ],"MyCandidatesForRemoval.txt",sep="\t", row.names=FALSE)
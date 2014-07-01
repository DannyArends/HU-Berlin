# Analysis of JAX known/pseudo gene data to create a custom Agilent array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#

setwd("E:/JAX/")

genedata <- read.table("MGI_Gene_Model_Coord.rpt", sep="\t", header=TRUE, row.names=1)
genesOnChr3 <- which(genedata[,"Ensembl.gene.chromosome"]==3)
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
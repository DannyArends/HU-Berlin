library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)

TFBSjasparALL <- query(MotifDb, "Mmusculus")

seqLogo(TFBSjasparALL[["Mmusculus-JASPAR_CORE-Hand1::Tcfe2a-MA0092.1"]])
seqLogo(TFBSjasparALL[["Mmusculus-cispb_1.02-M6097_1.02"]])
seqLogo(TFBSjasparALL[["Mmusculus-jolma2013-Elk3"]])
seqLogo(TFBSjasparALL[["Mmusculus-jolma2013-Elf5"]])
seqLogo(TFBSjasparALL[["Mmusculus-JASPAR_2014-GABPA-MA0062.2"]])
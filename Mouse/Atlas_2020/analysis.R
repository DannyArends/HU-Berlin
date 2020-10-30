setwd("D:/Edrive/Mouse/RNA/Atlas2020")

genotypes <- read.table("RAW/DN-2010_1861-Data/Analyse/1861_Genotypen.txt",sep="\t", header = TRUE,na.strings=c("NoCall", "", "---"), row.names=1)
annot <- read.table("RAW/MOUSEDIVm520650.na34.1.annot.csv/MOUSEDIVm520650.na34.1.annot.csv", sep = ",", header = TRUE, na.strings=c("NoCall", "", "---"), row.names=1)

datacols <- colnames(genotypes)[grep("chp.Call.Codes", colnames(genotypes), fixed=TRUE)]
gtsubset <- genotypes[, datacols]

datacols <- gsub("_20200702_.MOUSEDIVm520650..brlmm.p.chp.Call.Codes", "", gsub("X1861_", "", datacols))
datacols <- gsub("_20200707_.MOUSEDIVm520650..brlmm.p.chp.Call.Codes", "", datacols,fixed=TRUE)
datacols <- gsub("_20200703_.MOUSEDIVm520650..brlmm.p.chp.Call.Codes", "", datacols,fixed=TRUE)
datacols <- gsub(".", "-", datacols,fixed=TRUE)

colnames(gtsubset) <- datacols
gtsubset <- cbind(genotypes[,c(134:136)], annot[rownames(gtsubset), c("Allele.A", "Allele.B")], gtsubset)
write.table(gtsubset, "genotypesMouseDiversityArray.txt",sep = "\t", quote = FALSE, na = "")
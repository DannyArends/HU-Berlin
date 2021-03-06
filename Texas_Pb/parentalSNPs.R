setwd("D:/Edrive/Mouse/Texas_Pb/vcf")
for(f in list.files()[10:11]){
  outf <- paste0("../", gsub(".vcf", ".filtered.vcf", f))
  line.n  <- 1
  lengths <- c()
  Tfile <- file(f, "r")                                               # Get a pointer to the file
  cat("", file = outf)
  while(length((line = readLines(Tfile, n = 1) )) > 0){                                                   # Read a line, if available
    stsplit <- strsplit(line,"\t")[[1]]
    line.length <- length(stsplit)
    if(line.length == 11){
      cc0011 <- strsplit(stsplit[10], ":")[[1]][1]
      cc0017 <- strsplit(stsplit[11], ":")[[1]][1]
      if(cc0011 != "./." && cc0017 != "./." && cc0011 != cc0017){
        #cat("line:", line.n, "length:", line.length, " ",cc0011, " ",cc0017, "\n")                                # Number of words per line
        cat(line, "\n", file = outf, append=TRUE)
      }
    }
    lengths <- c(lengths, line.length)
    line.n <- line.n + 1                                                                                  # Increase out line number counter
  }
  close(Tfile)
}


setwd("D:/Edrive/Mouse/Texas_Pb/")
qtl2 <- read.csv("Texas_Pb_QTL2.filtered.vcf",sep="\t")
write.table(qtl2[,c(1:5, 10,11)], file = "Texas_Pb_QTL2.noinfo.vcf" , sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

setwd("D:/Edrive/Mouse/Texas_Pb/")
qtl2 <- read.csv("Texas_Pb_QTL3.filtered.vcf",sep="\t")
write.table(qtl2[,c(1:5, 10,11)], file = "Texas_Pb_QTL3.noinfo.vcf" , sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

setwd("D:/Edrive/Mouse/Texas_Pb/")
qtl2 <- read.csv("Texas_Pb_QTL4.filtered.vcf",sep="\t")
write.table(qtl2[,c(1:5, 10,11)], file = "Texas_Pb_QTL4.noinfo.vcf" , sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

setwd("D:/Edrive/Mouse/Texas_Pb/")
qtl2 <- read.csv("Texas_Pb_QTL5.filtered.vcf",sep="\t")
write.table(qtl2[,c(1:5, 10,11)], file = "Texas_Pb_QTL5.noinfo.vcf" , sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

setwd("D:/Edrive/Mouse/Texas_Pb/")
qtl2 <- read.csv("Texas_Pb_QTL6.filtered.vcf",sep="\t")
write.table(qtl2[,c(1:5, 10,11)], file = "Texas_Pb_QTL6.noinfo.vcf" , sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

setwd("D:/Edrive/Mouse/Texas_Pb/")
qtl2 <- read.csv("Texas_Pb_QTL9.filtered.vcf",sep="\t")
write.table(qtl2[,c(1:5, 10,11)], file = "Texas_Pb_QTL9.noinfo.vcf" , sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

setwd("D:/Edrive/Mouse/Texas_Pb/")
qtl2 <- read.csv("Texas_Pb_QTLI1-18.filtered.vcf",sep="\t")
write.table(qtl2[,c(1:5, 10,11)], file = "Texas_Pb_QTLI1-18.noinfo.vcf" , sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

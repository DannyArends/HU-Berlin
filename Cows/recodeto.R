
# Code from VCF coding to ref/alt coding
recodedVCF <- matrix(NA, nrow(seqVCF), 6 + length(10:ncol(seqVCF)), dimnames = list(c(), c("Chr", "Pos", "ChrPos", "ID", "Ref", "Alt", colnames(seqVCF)[10:ncol(seqVCF)])))
for(x in 1:nrow(seqVCF)){
  rowdata <- unlist(lapply(seqVCF[x, 10:ncol(seqVCF)], as.character))
  seqGT   <- unlist(lapply(strsplit(rowdata, ":"),"[",1))
  ref     <- seqVCF[x, "REF"]; alt <- seqVCF[x, "ALT"]
  hRef    <- paste0(ref,ref)
  hetero  <- paste0(sort(c(ref,alt)),collapse="")
  hAlt    <- paste0(alt,alt)
  seqGT[seqGT == "0/0"] <- hRef
  seqGT[seqGT == "0/1"] <- hetero
  seqGT[seqGT == "1/1"] <- hAlt
  seqGT[seqGT == "./."] <- NA
  recodedVCF[x,] <- c(Chr = seqVCF[x,"CHROM"], Pos = seqVCF[x,"POS"], ChrPos = paste(seqVCF[x,"CHROM"], seqVCF[x,"POS"], sep="_"), ID = seqVCF[x,"ID"], Ref = ref, Alt = alt, seqGT)
}

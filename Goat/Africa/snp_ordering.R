
getSNPorder <- function(data.all){
  snp.loc <- read.table("DNA/SihamAnalysis/FilteredLocationLWLT01.txt", sep = "\t", header = TRUE, row.names = 1)
  snp.loc <- cbind(snp.loc, pos = (snp.loc[, "Stop"] + snp.loc[, "Start"]) / 2)
  snp.loc <- snp.loc[sort(as.numeric(snp.loc[, "pos"]), index.return = TRUE)$ix,]
  snp.order <- c()
  for(x in 1:29){
    snp.order <- rbind(snp.order, snp.loc[which(snp.loc[, "chrN"] == x),])
  }
  snp.order <- snp.order[which(rownames(snp.order) %in% rownames(data.all)),]
  return(snp.order)
}

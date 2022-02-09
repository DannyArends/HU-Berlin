setwd("D:/Edrive/Mouse/Texas_Pb")

# Load the normalized counts
mdata <- read.table("NormalizedCounts.txt", header=TRUE, row.names=1, check.names=FALSE,na.strings=c("0","NA"))
mdata2 <- read.table("NormalizedCounts_4kidney.txt", header=TRUE, row.names=1, check.names=FALSE,na.strings=c("0","NA"))
mdata <- cbind(mdata, mdata2)

colnames(mdata) <- gsub(".sorted.bam","", colnames(mdata))
i <- which(colnames(mdata) == "16-CC011560-Liver")
colnames(mdata)[i] <- "16-CC011-560-Liver"
i <- which(colnames(mdata) == "28-CC011567-Kidney")
colnames(mdata)[i] <- "28-CC011-567-Kidney"
colnames(mdata) <- gsub("^.{0,3}", "", colnames(mdata))

samples.liver <- colnames(mdata)[grep("Liver", colnames(mdata))]
samples.kidney <- colnames(mdata)[grep("Kidney", colnames(mdata))]

library(preprocessCore)
ndata <- normalize.quantiles(as.matrix(log2(mdata + 1)))
rownames(ndata) <- rownames(mdata)
colnames(ndata) <- colnames(mdata)

#Liver: LDH, AST/LT, bilirubin and haptoglobin
#Kidney: hemoglobin

rows <- c(grep("Ldh", rownames(ndata)),grep("Ast", rownames(ndata)),grep("Hbb", rownames(ndata)),grep("Hp", rownames(ndata)))
ndata <- ndata[rows, ]

samples.liver <- colnames(mdata)[grep("Liver", colnames(mdata))]
samples.kidney <- colnames(mdata)[grep("Kidney", colnames(mdata))]

# Liver first
data.liver <- ndata[,samples.liver]
cc011.liver <- colnames(data.liver)[grep("CC011", colnames(data.liver))]
cc017.liver <- colnames(data.liver)[grep("CC017", colnames(data.liver))]

parental.liver <- ndata[,c(cc011.liver,cc017.liver)] # Take the 2 founders
parental.liver <- parental.liver[-which(apply(parental.liver,1, function(x){all(is.na(x))})),] # removed genes with all NA

DE.parental.liver <- apply(parental.liver,1,function(x){
  mean11 <- mean(x[1:length(cc011.liver)], na.rm=TRUE)
  sd11 <- sd(x[1:length(cc011.liver)], na.rm=TRUE)
  mean17 <- mean(x[(1+length(cc011.liver)): length(x)], na.rm=TRUE)
  sd17 <- sd(x[(1+length(cc011.liver)): length(x)], na.rm=TRUE)
  out <- tryCatch(t.test(x[1:length(cc011.liver)], x[(1+length(cc011.liver)): length(x)])$p.value, error=function(cond) {return(NA)})
  return(c(mean11, sd11, mean17, sd17, out))
})

liverSel <- cbind(round(t(DE.parental.liver)[,1:4],2), DE.parental.liver[5,])
colnames(liverSel) <- c("Mean(CC11)","SD(CC11)","Mean(CC17)", "SD(CC17)", "P")
liverSel <- liverSel[!is.na(liverSel[,5]),]

library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = "mgi_symbol", values = rownames(liverSel), mart = bio.mart)

rownames(res.biomart) <- res.biomart[, "mgi_symbol"]
liverSelAnn <- cbind(res.biomart[rownames(liverSel),], liverSel)
write.table(liverSelAnn, "Hemolysis.Liver.check.txt",sep="\t", quote=FALSE)
      

# kidney second
data.kidney <- ndata[,samples.kidney]
cc011.kidney <- colnames(data.kidney)[grep("CC011", colnames(data.kidney))]
cc017.kidney <- colnames(data.kidney)[grep("CC017", colnames(data.kidney))]

parental.kidney <- ndata[,c(cc011.kidney,cc017.kidney)] # Take the 2 founders
parental.kidney <- parental.kidney[-which(apply(parental.kidney,1, function(x){all(is.na(x))})),] # removed genes with all NA

DE.parental.kidney <- apply(parental.kidney,1,function(x){
  mean11 <- mean(x[1:length(cc011.kidney)], na.rm=TRUE)
  sd11 <- sd(x[1:length(cc011.kidney)], na.rm=TRUE)
  mean17 <- mean(x[(1+length(cc011.kidney)): length(x)], na.rm=TRUE)
  sd17 <- sd(x[(1+length(cc011.kidney)): length(x)], na.rm=TRUE)
  out <- tryCatch(t.test(x[1:length(cc011.kidney)], x[(1+length(cc011.kidney)): length(x)])$p.value, error=function(cond) {return(NA)})
  return(c(mean11, sd11, mean17, sd17, out))
})

kidneySel <- cbind(round(t(DE.parental.kidney)[,1:4],2), DE.parental.kidney[5,])
colnames(kidneySel) <- c("Mean(CC11)","SD(CC11)","Mean(CC17)", "SD(CC17)", "P")
kidneySel <- kidneySel[!is.na(kidneySel[,5]),]

library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = "mgi_symbol", values = rownames(kidneySel), mart = bio.mart)

rownames(res.biomart) <- res.biomart[, "mgi_symbol"]
kidneySelAnn <- cbind(res.biomart[rownames(kidneySel),], kidneySel)
write.table(kidneySelAnn, "Hemolysis.Kidney.check.txt",sep="\t", quote=FALSE)

      
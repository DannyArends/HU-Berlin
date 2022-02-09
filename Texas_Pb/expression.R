setwd("D:/Edrive/Mouse/Texas_Pb")

regions <- c("1:124548738:184926264", 
             "1:135464050:157672910", 
             "3:85146443:147514132", 
             "4:58831312:88222686", 
             "6:76940963:114257130", 
             "7:62305639:81923457", 
             "8:31799796:104944836",
             "1:151111612:157672910",
             "18:81259093:88823124")

# Load the genes in the QTL 95% CI
res <- vector("list", length(regions))
i <- 1
for(r in regions){
  res[[i]] <- read.table(file=paste0("genes_", gsub(":", "-",r), ".txt"), sep="\t", header=TRUE)
  cat("Done",r,"\n")
  i <- i + 1
}

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

heatmap(cor(ndata, use="pair"))

# Do the liver first
data.liver <- ndata[,samples.liver]
cc011.liver <- colnames(data.liver)[grep("CC011", colnames(data.liver))]
cc017.liver <- colnames(data.liver)[grep("CC017", colnames(data.liver))]

parental.liver <- ndata[,c(cc011.liver,cc017.liver)] # Take the 2 founders
parental.liver <- parental.liver[-which(apply(parental.liver,1, function(x){all(is.na(x))})),] # removed genes with all NA

DE.parental.liver <- apply(parental.liver,1,function(x){
  out <- tryCatch(t.test(x[1:length(cc011.liver)], x[(1+length(cc011.liver)): length(x)])$p.value, error=function(cond) {return(NA)})
  return(out)
})

plot(-log10(DE.parental.liver))
abline(h = -log10(0.01 / length(DE.parental.liver)))

significant.liver <- names(which(-log10(DE.parental.liver) > -log10(0.05 / length(DE.parental.liver))))

# Do the kidney second
data.kidney <- ndata[,samples.kidney]
cc011.kidney <- colnames(data.kidney)[grep("CC011", colnames(data.kidney))]
cc017.kidney <- colnames(data.kidney)[grep("CC017", colnames(data.kidney))]

parental.kidney <- ndata[,c(cc011.kidney,cc017.kidney)] # Take the 2 founders
parental.kidney <- parental.kidney[-which(apply(parental.kidney,1, function(x){all(is.na(x))})),] # removed genes with all NA

DE.parental.kidney <- apply(parental.kidney,1,function(x){
  out <- tryCatch(t.test(x[1:length(cc011.kidney)], x[(1+length(cc011.kidney)): length(x)])$p.value, error=function(cond) {return(NA)})
  return(out)
})

plot(-log10(DE.parental.kidney))
abline(h = -log10(0.01 / length(DE.parental.kidney)))

significant.kidney <- names(which(-log10(DE.parental.kidney) > -log10(0.05 / length(DE.parental.kidney))))
#unknown <- c(grep("^Gm", names(significant)), grep("^LOC", names(significant)), grep("Rik$", names(significant)))
#significant <- names(significant[-unknown])

library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = "mgi_symbol", values = c(significant.liver, significant.kidney), mart = bio.mart)

inRegionDiffEx <- NULL
for(r in 1:length(regions)){
  subS <- res.biomart[which(res.biomart[, "ensembl_gene_id"] %in% res[[r]][, "ensembl_gene_id"]),]
  Ecc11.liver <- parental.liver[subS[,"external_gene_name"],1:length(cc011.liver)]
  Ecc11.kidney <- parental.kidney[subS[,"external_gene_name"],1:length(cc011.kidney)]

  Ecc17.liver <- parental.liver[subS[,"external_gene_name"],(1+length(cc011.liver)): ncol(parental.liver)]
  Ecc17.kidney <- parental.kidney[subS[,"external_gene_name"],(1+length(cc011.kidney)): ncol(parental.kidney)]
  if(nrow(subS) > 1){
    cat("1+ genes\n")
    mcc11.liver <- round(as.numeric(apply(Ecc11.liver,1,mean,na.rm=TRUE)),2)
    mcc11.kidney <- round(as.numeric(apply(Ecc11.kidney,1,mean,na.rm=TRUE)),2)
    sdcc11.liver <- round(as.numeric(apply(Ecc11.liver,1,sd,na.rm=TRUE)),2)
    sdcc11.kidney <- round(as.numeric(apply(Ecc11.kidney,1,sd,na.rm=TRUE)),2)

    mcc17.liver <- round(as.numeric(apply(Ecc17.liver,1,mean,na.rm=TRUE)),2)
    mcc17.kidney <- round(as.numeric(apply(Ecc17.kidney,1,mean,na.rm=TRUE)),2)
    sdcc17.liver <- round(as.numeric(apply(Ecc17.liver,1,sd,na.rm=TRUE)),2)
    sdcc17.kidney <- round(as.numeric(apply(Ecc17.kidney,1,sd,na.rm=TRUE)),2)
    inRegionDiffEx <- rbind(inRegionDiffEx, cbind(region = regions[r], subS, 
                                                  Liver_CC011 = paste0(mcc11.liver," (", sdcc11.liver, ")"), Liver_CC017 = paste0(mcc17.liver," (", sdcc17.liver, ")"),
                                                  Liver_P = DE.parental.liver[subS[,"external_gene_name"]],
                                                  Kidney_CC011 = paste0(mcc11.kidney," (", sdcc11.kidney, ")"), Kidney_CC017 = paste0(mcc17.kidney," (", sdcc17.kidney, ")"),
                                                  Kidney_P = DE.parental.kidney[subS[,"external_gene_name"]]
                                                  ))
  }
  if(nrow(subS) == 1){
    cat("1 gene\n")
    mcc11.liver <- round(mean(Ecc11.liver,na.rm=TRUE),2)
    mcc11.kidney <- round(mean(Ecc11.kidney,na.rm=TRUE),2)
    sdcc11.liver <- round(sd(Ecc11.liver,na.rm=TRUE),2)
    sdcc11.kidney <- round(sd(Ecc11.kidney,na.rm=TRUE),2)
    mcc17.liver <- round(mean(Ecc17.liver,na.rm=TRUE),2)
    mcc17.kidney <- round(mean(Ecc17.kidney,na.rm=TRUE),2)
    sdcc17.liver <- round(sd(Ecc17.liver,na.rm=TRUE),2)
    sdcc17.kidney <- round(sd(Ecc17.kidney,na.rm=TRUE),2)
    inRegionDiffEx <- rbind(inRegionDiffEx, unlist(c(region = regions[r], subS, 
                                                  Liver_CC011 = paste0(mcc11.liver," (", sdcc11.liver, ")"), Liver_CC017 = paste0(mcc17.liver," (", sdcc17.liver, ")"), Liver_P = DE.parental.liver[subS[,"external_gene_name"]],
                                                  Kidney_CC011 = paste0(mcc11.kidney," (", sdcc11.kidney, ")"), Kidney_CC017 = paste0(mcc17.kidney," (", sdcc17.kidney, ")"), Kidney_P = DE.parental.kidney[subS[,"external_gene_name"]]
                                                  )))
  }
}
rownames(inRegionDiffEx) <- 1:nrow(inRegionDiffEx)
write.table(inRegionDiffEx, "diffExpression.parental.txt", sep = "\t", quote=FALSE, row.names=FALSE)


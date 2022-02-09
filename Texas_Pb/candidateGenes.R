setwd("D:/Edrive/Mouse/Texas_Pb")

library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="http://nov2020.archive.ensembl.org")

regions <- c("1:124548738:184926264", 
             "1:135464050:157672910", 
             "3:85146443:147514132", 
             "4:58831312:88222686", 
             "6:76940963:114257130", 
             "7:62305639:81923457", 
             "7:62305639:97260868", 
             "7:19482853:119720011", 
             "8:31799796:104944836",
             "1:151111612:157672910",
             "18:81259093:88823124")

for(r in regions[10:11]){
  res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("chromosomal_region", "biotype"), values = list(r, "protein_coding"), mart = bio.mart)
  write.table(res.biomart, file=paste0("genes_", gsub(":", "-",r), ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
  cat("Done",r,"\n")
}

# VEP predictions and CC founder genotypes
setwd("D:/Edrive/Mouse/Texas_Pb")
library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="http://nov2020.archive.ensembl.org")
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("1:135464050:157672910", "protein_coding"), mart = bio.mart)


QTL2.founder <- read.table("CCfounders_QTL2.snps-filtered.vcf", sep = "\t")
colnames(QTL2.founder) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","A_J","C57BL_6NJ","129S1_SvImJ","NOD_ShiLtJ","NZO_HlLtJ","CAST_EiJ","PWK_PhJ","WSB_EiJ")
#PWK allele should be the one with the effect
QTL2.founder[, "A_J"] <- unlist(lapply(strsplit(QTL2.founder[,"A_J"], ":"), "[", 1))
QTL2.founder[, "C57BL_6NJ"] <- unlist(lapply(strsplit(QTL2.founder[,"C57BL_6NJ"], ":"), "[", 1))
QTL2.founder[, "129S1_SvImJ"] <- unlist(lapply(strsplit(QTL2.founder[,"129S1_SvImJ"], ":"), "[", 1))
QTL2.founder[, "NOD_ShiLtJ"] <- unlist(lapply(strsplit(QTL2.founder[,"NOD_ShiLtJ"], ":"), "[", 1))
QTL2.founder[, "NZO_HlLtJ"] <- unlist(lapply(strsplit(QTL2.founder[,"NZO_HlLtJ"], ":"), "[", 1))
QTL2.founder[, "CAST_EiJ"] <- unlist(lapply(strsplit(QTL2.founder[,"CAST_EiJ"], ":"), "[", 1))
QTL2.founder[, "PWK_PhJ"] <- unlist(lapply(strsplit(QTL2.founder[,"PWK_PhJ"], ":"), "[", 1))
QTL2.founder[, "WSB_EiJ"] <- unlist(lapply(strsplit(QTL2.founder[,"WSB_EiJ"], ":"), "[", 1))
QTL2.founder <- QTL2.founder[, - c(3,6,7,8,9)]
keep <- c()
for(x in 1:nrow(QTL2.founder)){
  if(!QTL2.founder[x, "PWK_PhJ"] %in% QTL2.founder[x, c("C57BL_6NJ", "129S1_SvImJ", "NOD_ShiLtJ","NZO_HlLtJ","C57BL_6NJ")]){
    keep <- c(keep, x)
  }
}
QTL2.founder <- QTL2.founder[keep,]

#QTL2.parentals <- read.table("Texas_Pb_QTL2.filtered.vcf", sep = "\t")
#colnames(QTL2.parentals) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "CC011","CC017")

QTL2.vep <- read.table("Texas_Pb_QTL2.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
QTL2.vep <- QTL2.vep[which(QTL2.vep[, "V5"] %in% c("HIGH", "MODERATE")),]
QTL2.vep <- QTL2.vep[-which(duplicated(QTL2.vep[, "V2"])),]
QTL2.summ <- table(QTL2.vep[,"V6"])



PWK.unique <- paste0(QTL2.founder[,"CHROM"], ":", QTL2.founder[,"POS"], "-", QTL2.founder[,"POS"])
QTL2.vep <- QTL2.vep[which(QTL2.vep[,2] %in% PWK.unique),]
QTL2.pwk.summ <- table(QTL2.vep[,"V6"])
 
res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA, "SNP from PWK_PhJ" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTL2.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTL2.summ[ii]
  }
  ii <- which(names(QTL2.pwk.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "SNP from PWK_PhJ"] = 0
  }else{
    res.biomart[x, "SNP from PWK_PhJ"] = QTL2.pwk.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTL2.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)



# VEP of QTL 1
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("1:124548738:184926264", "protein_coding"), mart = bio.mart)

QTL1.vep <- read.table("Texas_Pb_QTL1.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
imp <- unlist(lapply(strsplit(QTL1.vep[,14], ";"),"[",1))
QTL1.vep <- QTL1.vep[which(imp %in% c("IMPACT=HIGH", "IMPACT=MODERATE")),]
QTL1.vep <- QTL1.vep[-which(duplicated(QTL1.vep[, "V2"])),]
QTL1.summ <- table(QTL1.vep[,"V4"])

res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTL1.summ) == res.biomart[x, "ensembl_gene_id"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTL1.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTL1.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)





# VEP predictions and CC founder genotypes
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("3:85146443:147514132", "protein_coding"), mart = bio.mart)

QTL3.vep <- read.table("Texas_Pb_QTL3.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
QTL3.vep <- QTL3.vep[which(QTL3.vep[, "V5"] %in% c("HIGH", "MODERATE")),]
QTL3.vep <- QTL3.vep[-which(duplicated(QTL3.vep[, "V2"])),]
QTL3.summ <- table(QTL3.vep[,"V6"])

res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTL3.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTL3.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTL3.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)


# VEP of QTL 4 
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("4:58831312:88222686", "protein_coding"), mart = bio.mart)

QTL4.vep <- read.table("Texas_Pb_QTL4.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
QTL4.vep <- QTL4.vep[which(QTL4.vep[, "V5"] %in% c("HIGH", "MODERATE")),]
QTL4.vep <- QTL4.vep[-which(duplicated(QTL4.vep[, "V2"])),]
QTL4.summ <- table(QTL4.vep[,"V6"])

res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTL4.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTL4.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTL4.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)


# VEP of QTL 5
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("6:76940963:114257130", "protein_coding"), mart = bio.mart)

QTL5.vep <- read.table("Texas_Pb_QTL5.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
QTL5.vep <- QTL5.vep[which(QTL5.vep[, "V5"] %in% c("HIGH", "MODERATE")),]
QTL5.vep <- QTL5.vep[-which(duplicated(QTL5.vep[, "V2"])),]
QTL5.summ <- table(QTL5.vep[,"V6"])

res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTL5.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTL5.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTL5.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)


# VEP of QTL 6
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("7:62305639:81923457", "protein_coding"), mart = bio.mart)

QTL6.vep <- read.table("Texas_Pb_QTL6.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
QTL6.vep <- QTL6.vep[which(QTL6.vep[, "V5"] %in% c("HIGH", "MODERATE")),]
QTL6.vep <- QTL6.vep[-which(duplicated(QTL6.vep[, "V2"])),]
QTL6.summ <- table(QTL6.vep[,"V6"])

res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTL6.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTL6.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTL6.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)


# VEP of QTL 9
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("8:31799796:104944836", "protein_coding"), mart = bio.mart)

QTL9.vep <- read.table("Texas_Pb_QTL9.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
imp <- unlist(lapply(strsplit(QTL9.vep[,14], ";"),"[",1))
QTL9.vep <- QTL9.vep[which(imp %in% c("IMPACT=HIGH", "IMPACT=MODERATE")),]
QTL9.vep <- QTL9.vep[-which(duplicated(QTL9.vep[, "V2"])),]
QTL9.summ <- table(QTL9.vep[,"V4"])

res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTL9.summ) == res.biomart[x, "ensembl_gene_id"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTL9.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTL9.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)



# VEP predictions and CC founder genotypes chromosome 1 of the 1-18 interaction
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("1:151111612:157672910", "protein_coding"), mart = bio.mart)

QTLI1.vep <- read.table("Texas_Pb_QTLI1-1.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
QTLI1.vep <- QTLI1.vep [c(grep("HIGH", QTLI1.vep[,"V14"]), grep("MODERATE", QTLI1.vep[,"V14"])), ]
QTLI1.vep <- QTLI1.vep[-which(duplicated(QTLI1.vep[, "V2"])),]

seppVals <- strsplit(QTLI1.vep[,"V14"], ";")
symbolP <- unlist(lapply(seppVals, function(x){grep("SYMBOL=",x, fixed=TRUE)}))
gsymbols <- c()
for(x in 1:length(seppVals)){
  gsymbols <- c(gsymbols, gsub("SYMBOL=", "", seppVals[[x]][symbolP[x]]))
}

QTLI1.summ <- table(gsymbols)

res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTLI1.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTLI1.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTLI1-1.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)



# VEP predictions and CC founder genotypes chromosome 18 of the 1-18 interaction
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), values = list("18:81259093:88823124", "protein_coding"), mart = bio.mart)

QTLI1.vep <- read.table("Texas_Pb_QTLI1-18.vep.txt", sep = "\t")
# Take high and moderate impact SNPs
QTLI1.vep <- QTLI1.vep [c(grep("HIGH", QTLI1.vep[,"V14"]), grep("MODERATE", QTLI1.vep[,"V14"])), ]
QTLI1.vep <- QTLI1.vep[-which(duplicated(QTLI1.vep[, "V2"])),]

seppVals <- strsplit(QTLI1.vep[,"V14"], ";")
symbolP <- unlist(lapply(seppVals, function(x){grep("SYMBOL=",x, fixed=TRUE)}))
gsymbols <- c()
for(x in 1:length(seppVals)){
  gsymbols <- c(gsymbols, gsub("SYMBOL=", "", seppVals[[x]][symbolP[x]]))
}

QTLI1.summ <- table(gsymbols)


res.biomart <- cbind(res.biomart, "Impactful SNPs" = NA)
for(x in 1:nrow(res.biomart)){
  ii <- which(names(QTLI1.summ) == res.biomart[x, "mgi_symbol"])
  if(length(ii) == 0){
    res.biomart[x, "Impactful SNPs"] = 0
  }else{
    res.biomart[x, "Impactful SNPs"] = QTLI1.summ[ii]
  }
}

write.table(res.biomart, file=paste0("candidateGenes/QTLI1-18.annot.vep.txt"), sep="\t", quote=FALSE, row.names=FALSE)


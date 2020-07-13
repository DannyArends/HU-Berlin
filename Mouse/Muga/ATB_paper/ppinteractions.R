setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
results <- read.table("ChiSqScores_Incompatible.txt", sep="\t", check.names=FALSE)
rr <- apply(results,1,as.numeric)
rownames(rr) <- colnames(rr)
results <- rr

# Load the Allele Transmission Biased regions
regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)

allRegions <- rbind(cbind(regionsMat, origin = "MAT"), cbind(regionsPat, origin = "PAT"))
allRegions <- allRegions[with(allRegions, order(Chr, Start)), ]

regionnames <- paste0(allRegions[,"Chr"],":", 
          round(allRegions[,"Start"] / 1000000, 0),"-",
          round(allRegions[,"Stop"] / 1000000, 0), " ", 
          allRegions[,"origin"])

rownames(allRegions) <- regionnames

allRegions["8:62-63 PAT", "Prefered.Allele"] <- "B6N"
allRegions["9:88-91 MAT", "Prefered.Allele"] <- "B6N"
allRegions["11:13-13 MAT", "Prefered.Allele"] <- "BFMI"


# Generate an overview plot
nTests <- (ncol(results) * ncol(results)) / 2
LODscores <- -log10(pchisq(results, 1, lower.tail=FALSE))
threshold <- -log10(0.05 / nTests)

significant <-  names(which(apply(LODscores,1,max,na.rm=TRUE) > threshold))

allRegions <- allRegions[significant, ]
results <- results[significant,significant]

bfmipref <- which(allRegions[,"Prefered.Allele"] == "BFMI")
b6npref <- which(allRegions[,"Prefered.Allele"] == "B6N")

pat_bfmi_snp <- read.table(file="Analysis/PAT_BFMI_NSyn_Filtered.txt", sep = "\t", header=TRUE, colClasses="character")
mat_bfmi_snp <- read.table(file="Analysis/MAT_BFMI_NSyn_Filtered.txt", sep = "\t", header=TRUE, colClasses="character")
pat_b6_snp <- read.table(file="Analysis/PAT_B6N_NSyn_Filtered.txt", sep = "\t", header=TRUE, colClasses="character")
mat_b6_snp <- read.table(file="Analysis/MAT_B6N_NSyn_Filtered.txt", sep = "\t", header=TRUE, colClasses="character")

all_snps <- rbind(pat_bfmi_snp, mat_bfmi_snp, pat_b6_snp, mat_b6_snp)[, 1:6]
all_snps <- all_snps[-which(duplicated(all_snps)),]

length(unique(c(as.character(pat_bfmi_snp[,1]), as.character(mat_bfmi_snp[,1]), as.character(pat_b6_snp[,1]), as.character(mat_b6_snp[,1]))),sep="\n")

cat(unique(c(as.character(pat_bfmi_snp[,5]), as.character(mat_bfmi_snp[,5]), as.character(pat_b6_snp[,5]), as.character(mat_b6_snp[,5]))),sep="\n")

known_interactions <- read.csv("InnateDB_genes.txt", sep = "\t", colClasses="character")
known_interactions <- known_interactions[which(known_interactions[, "ensembl"] %in% all_snps[,1]),]
dim(known_interactions)
onames <- c()
for(g in known_interactions[, "QUERY.XREF"]){
  onames <- c(onames, all_snps[which(all_snps[,1] == g),5])
}
known_interactions <- cbind(oname = onames, known_interactions)
known_interactions <- known_interactions[, colnames(known_interactions) %in% c("oname", "name", "nrIntxsValidated")]
dim(known_interactions)[1] /2

all_snps[which(all_snps[,5] == "Asph"),]
all_snps[which(all_snps[,5] == "ApoB"),]
all_snps[which(all_snps[,5] == "Hs1bp3"),]
all_snps[which(all_snps[,5] == "Akap2"),]

allgenes <- c()
for(x in rownames(results)){
  cat("------------------------\n")
  splitted <- strsplit(x, ":")[[1]]
  chrO <- splitted[1]
  posO <- strsplit(splitted[2], " ")[[1]][1]
  posOs <- (as.numeric(strsplit(posO, "-")[[1]][1])-1) * 1000000
  posOe <- (as.numeric(strsplit(posO, "-")[[1]][2])+1) * 1000000
  inR <- which(all_snps[, "chromosome_name"] == chrO & all_snps[, "start_position"] > posOs & all_snps[, "end_position"] < posOe)
  if(length(inR) > 0) {
    genes <- as.character(all_snps[inR, "ensembl_gene_id"])
    cat(x, chrO, posOs, posOe, " ", length(genes), "\n")
    allgenes <- c(allgenes, genes)
  }
}

allgenes <- unique(allgenes)


library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

res.biomart <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("ensembl_gene_id"), values = all_snps[,1], mart = bio.mart)

mdata <- read.table("10090.protein.links.v11.0.txt", header = TRUE)

mdata[,1] <- gsub("10090.", "", mdata[,1], fixed=TRUE)
mdata[,2] <- gsub("10090.", "", mdata[,2], fixed=TRUE)

mdata <- mdata[which(mdata[,1] %in% res.biomart[, "ensembl_peptide_id"]),]
mdata <- mdata[which(mdata[,2] %in% res.biomart[, "ensembl_peptide_id"]),]

mdatafull <- c()
for(x in 1:nrow(mdata)){
  p1gene <- res.biomart[which(res.biomart[,2] == mdata[x,1]),c("ensembl_gene_id", "mgi_symbol")]
  p2gene <- res.biomart[which(res.biomart[,2] == mdata[x,2]),c("ensembl_gene_id", "mgi_symbol")]
  mdatafull <- rbind(mdatafull, c(mdata[x,1], p1gene,mdata[x,2],p2gene,mdata[x,3]))
}
colnames(mdatafull) <- c("protein1", "gene1", "symbol1", "protein2", "gene2", "symbol2", "score")



res.biomart.region <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("ensembl_gene_id"), values = allgenes, mart = bio.mart)

mdata <- read.table("10090.protein.links.v11.0.txt", header = TRUE)

mdata[,1] <- gsub("10090.", "", mdata[,1], fixed=TRUE)
mdata[,2] <- gsub("10090.", "", mdata[,2], fixed=TRUE)

mdata <- mdata[which(mdata[,1] %in% res.biomart.region[, "ensembl_peptide_id"]),]
mdata <- mdata[which(mdata[,2] %in% res.biomart.region[, "ensembl_peptide_id"]),]

mdatafull.region <- c()
for(x in 1:nrow(mdata)){
  p1gene <- res.biomart.region[which(res.biomart.region[,2] == mdata[x,1]),c("ensembl_gene_id", "mgi_symbol")]
  p2gene <- res.biomart.region[which(res.biomart.region[,2] == mdata[x,2]),c("ensembl_gene_id", "mgi_symbol")]
  mdatafull.region <- rbind(mdatafull.region, c(mdata[x,1], p1gene,mdata[x,2],p2gene,mdata[x,3]))
}
colnames(mdatafull.region) <- c("protein1", "gene1", "symbol1", "protein2", "gene2", "symbol2", "score")



allgenes <- c()
for(x in rownames(results)){
  cat("------------------------\n")
  fp <- paste0("Analysis/", allRegions[x, "origin"], "_", allRegions[x, "Prefered.Allele"], "_NSyn_Filtered.txt")
  splitted <- strsplit(x, ":")[[1]]
  chrO <- splitted[1]
  posO <- strsplit(splitted[2], " ")[[1]][1]
  posOs <- (as.numeric(strsplit(posO, "-")[[1]][1])-1) * 1000000
  posOe <- (as.numeric(strsplit(posO, "-")[[1]][2])+1) * 1000000
  inR <- which(all_snps[, "chromosome_name"] == chrO & all_snps[, "start_position"] > posOs & all_snps[, "end_position"] < posOe)
  if(length(inR) > 0) {
    genes <- as.character(all_snps[inR, "mgi_symbol"])
    cat(x, " ", length(genes), "\n")
    for(g1 in genes){
      allgenes <- c(allgenes, g1)
      for(y in names(which(results[x,] > threshold))){
        fp2 <- paste0("Analysis/", allRegions[y, "origin"], "_", allRegions[y, "Prefered.Allele"], "_NSyn_Filtered.txt")
        splitted <- strsplit(y, ":")[[1]]
        chrR <- splitted[1]
        posR <- strsplit(splitted[2], " ")[[1]][1]
        posRs <- (as.numeric(strsplit(posR, "-")[[1]][1])-1) * 1000000
        posRe <- (as.numeric(strsplit(posR, "-")[[1]][2])+1) * 1000000
        inR2 <- which(all_snps[, "chromosome_name"] == chrR & all_snps[, "start_position"] > posRs & all_snps[, "end_position"] < posRe)
        if(length(inR2) > 0) {
          genes2 <- as.character(all_snps[inR2, "mgi_symbol"])
          for(g2 in genes2){
            allgenes <- c(allgenes, g2)
            ii <- which(known_interactions[,1] == g1 & known_interactions[,2] == g2)
            if(length(ii) > 0) cat(x, g1, y, g2, known_interactions[ii,3], "\n")
          }
        }
      }
    }
  }
}





#http://mips.helmholtz-muenchen.de/proj/ppi/
library("XML")
ppi_xml <- xmlParse(file = "allppis.xml")
rootnode <- xmlRoot(ppi_xml)
xmllist <- xmlToList(rootnode)



#  sp: Swissprot and TREMBL
#  gb: Genbank
#  pir: PIR

ppi_table <- c()
for(x in 1:length(xmllist[[1]]$interactionList)){
  ref <- paste0(xmllist[[1]]$interactionList[[x]]$experimentList$experimentDescription$bibref$xref$primaryRef,collapse=":")
  expi <- xmllist[[1]]$interactionList[[x]]$experimentList$experimentDescription$interactionDetection$names$shortLabel
  A <- xmllist[[1]]$interactionList[[x]]$participantList[[1]]$proteinInteractor$names$fullName
  AID <- paste0(xmllist[[1]]$interactionList[[x]]$participantList[[1]]$proteinInteractor$xref$primaryRef, collapse=":")
  B <- xmllist[[1]]$interactionList[[x]]$participantList[[2]]$proteinInteractor$names$fullName
  BID <- paste0(xmllist[[1]]$interactionList[[x]]$participantList[[2]]$proteinInteractor$xref$primaryRef, collapse=":")
  ppi_table <- rbind(ppi_table, c(A,AID,B, BID, ref, expi))
}
colnames(ppi_table) <- c("proteinA", "DBIDA", "proteinB", "DBIDB", "Publication", "Type")

library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("uniprotsptrembl"), values = list(r, "protein_coding"), mart = bio.mart)

getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), values = "U47050", mart = ensembl)
                     
for(x in rownames(results)){
  genes <- NULL
  cat("------------------------\n")
  fp <- paste0("Analysis/", allRegions[x, "origin"], "_", allRegions[x, "Prefered.Allele"], "_NSyn_Filtered.txt")
  geneTableO <- read.table(file = fp, sep = "\t", header=TRUE)
  splitted <- strsplit(x, ":")[[1]]
  chrO <- splitted[1]
  posO <- strsplit(splitted[2], " ")[[1]][1]
  posOs <- (as.numeric(strsplit(posO, "-")[[1]][1])-1) * 1000000
  posOe <- (as.numeric(strsplit(posO, "-")[[1]][2])+1) * 1000000
  inR <- which(geneTableO[, "chromosome_name"] == chrO & geneTableO[, "start_position"] > posOs & geneTableO[, "end_position"] < posOe)
  if(length(inR) > 0) {
    genes <- as.character(geneTableO[inR, "mgi_symbol"])
    for(g in genes){
      cat(x, g, " ", "*\n")
    }
    for(y in names(which(results[x,] > threshold))){
      fp2 <- paste0("Analysis/", allRegions[y, "origin"], "_", allRegions[y, "Prefered.Allele"], "_NSyn_Filtered.txt")
      geneTableR <- read.table(file = fp, sep = "\t", header=TRUE)
      splitted <- strsplit(y, ":")[[1]]
      chrR <- splitted[1]
      posR <- strsplit(splitted[2], " ")[[1]][1]
      posRs <- (as.numeric(strsplit(posR, "-")[[1]][1])-1) * 1000000
      posRe <- (as.numeric(strsplit(posR, "-")[[1]][2])+1) * 1000000
      inR2 <- which(geneTableO[, "chromosome_name"] == chrR & geneTableO[, "start_position"] > posRs & geneTableO[, "end_position"] < posRe)
      if(length(inR2) > 0) {
        genes2 <- as.character(geneTableR[inR2, "mgi_symbol"])
        for(g in genes2){
          cat(y, g, " ", "+\n")
        }
      }
    }
  }
  cat("------------------------\n")
}


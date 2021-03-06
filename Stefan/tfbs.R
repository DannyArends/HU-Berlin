# TFBS
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015
library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
#RPKM <- read.table(file="Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep="\t")
#SUBSET <- RPKM[which(RPKM[,"Ratio_F1_Scale"] > 1.2 | RPKM[,"Ratio_F1_Scale"] < -1.2),]
#SUBSET <- SUBSET[which(SUBSET[,"tTest_F1"] < 0.1),]
#SUBSET <- SUBSET[which(SUBSET[,"Mean.BFMI860.12xB6N.L"] > 0.5 | SUBSET[,"Mean.B6NxBFMI860.12.L"] > 0.5),]
#genes <- SUBSET[,"ensembl_gene_id"]
#genes <- as.character(unlist(read.table("S2_ENS.txt")))
genes <- c("ENSMUSG00000038521","ENSMUSG00000101795","ENSMUSG00000100774","ENSMUSG00000032235","ENSMUSG00000031604","ENSMUSG00000032786","ENSMUSG00000101514","ENSMUSG00000059824","ENSMUSG00000028680","ENSMUSG00000032078","ENSMUSG00000031431","ENSMUSG00000101517","ENSMUSG00000076617",
"ENSMUSG00000025203","ENSMUSG00000028712","ENSMUSG00000021228","ENSMUSG00000079343","ENSMUSG00000058230","ENSMUSG00000003464","ENSMUSG00000079507","ENSMUSG00000066477","ENSMUSG00000098623","ENSMUSG00000024222","ENSMUSG00000006711","ENSMUSG00000032376","ENSMUSG00000060068","ENSMUSG00000019132",
"ENSMUSG00000050621","ENSMUSG00000052632","ENSMUSG00000005069","ENSMUSG00000028445","ENSMUSG00000062901","ENSMUSG00000025369","ENSMUSG00000074115","ENSMUSG00000028494","ENSMUSG00000029254","ENSMUSG00000044254","ENSMUSG00000083365","ENSMUSG00000041828","ENSMUSG00000025202","ENSMUSG00000030555",
"ENSMUSG00000022814","ENSMUSG00000027109","ENSMUSG00000024870","ENSMUSG00000034645","ENSMUSG00000069267","ENSMUSG00000025791","ENSMUSG00000041935","ENSMUSG00000042506","ENSMUSG00000036478","ENSMUSG00000034135","ENSMUSG00000031980","ENSMUSG00000040351","ENSMUSG00000035372","ENSMUSG00000036305",
"ENSMUSG00000049232","ENSMUSG00000064181","ENSMUSG00000031196","ENSMUSG00000040562","ENSMUSG00000072115","ENSMUSG00000032120","ENSMUSG00000018340","ENSMUSG00000038641","ENSMUSG00000028307","ENSMUSG00000024371","ENSMUSG00000095597","ENSMUSG00000039156","ENSMUSG00000022550")
genes <- genes[-which(!genes %in% names(transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene")))]
chromosomal.loc <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene") [genes]
promoter.gene <- getPromoterSeq(chromosomal.loc, Mmusculus, upstream=2000, downstream=1000)
TFBSjasparALL <- query(MotifDb, "Mmusculus")

# nur die aus estrogen (oder auswahl)
estrogen <- c("Mmusculus-jolma2013-Esrra-2", "Mmusculus-UniPROBE-Esrra.UP00079", "Mmusculus-JASPAR_CORE-Esrrb-MA0141.1", "Mmusculus-JASPAR_2014-Esrrb-MA0141.2","Mmusculus-jolma2013-Ar")
estrogen <- c("Mmusculus-jolma2013-Ar")

TFBSjasparALL <- TFBSjasparALL[estrogen]
output <- matrix(NA, length(names(TFBSjasparALL)),length(genes), dimnames=list(names(TFBSjasparALL), genes))

### just estrogen Esrr
for(tfbs in 1:length(TFBSjasparALL)){
seqLogo(TFBSjasparALL[[tfbs]])
TFBSjaspar <- round(100 * TFBSjasparALL[[tfbs]])
for(gene in names(promoter.gene)){
hits <- matchPWM(TFBSjaspar, unlist(promoter.gene[gene])[[1]], "90%")
output[names(TFBSjasparALL)[tfbs], gene] <- length(hits)
}
}

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") # Biomart for mouse genes
fullnames <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), filters = c("ensembl_gene_id"), values = colnames(output), mart = bio.mart)
output2 <- output[,match(fullnames[,"ensembl_gene_id"], colnames(output))]
write.table(output, file="TFBSmatches_PaperGenes_Ar.90prz.txt",sep="\t")
write.table(rbind(fullnames[,"mgi_symbol"], output2), file="TFBSmatches_PaperGenes_Ar.90prz.txt",sep="\t",quote=FALSE)



#### table 3 hoghest p-values + Cyp
hi <- c("Slc16a7","Ehhadh","Acox1","Ccny","Xdh","Pla2g12b","Trmt2b","Cul4a","Abcd3","Aplp2","Bloc1s5","Il1r1","Rapgef4","Rbbp4","Rarres2","Cela1","Ssr4","Aco1","Mtus1","Rpa1","Rpl32","Otud7b","Nlrp12","Rps4x","C1s","Rassf6","Hbp1","Mgll","Cpped1","Coq10b","Gm11512","Gm13340","Aldh3a2","Ivd","Pts","Gm5835","Vps4b", "Cyp4a10","Cyp4a14")
hi <- c("Peg3","Zrsr1","H13","Agpat9","Arhgef40","Gbp11","Gm4735","Gm8893","Lbp","Serpina11","Tpm3-rs7","Ttc1","C4b","Ceacam1","Ces1g","Cfhr1","Clns1a","Cyp2c37","Cyp2c68","Cyp4a14","Dnm1l","Dusp12","Gm13775","Gm2382","H2-K1","Hist1h4f","Iah1","Mug2","Nnmt","Osgep","Rbm12b1","Rpl29","Serpina12","Serpina1a","Serpina1d","Sfr1","Slc22a28","Tma7","Vnn1")
hi <- RPKM[which(RPKM$mgi_symbol %in% hi),"ensembl_gene_id"]
genes <- hi

chromosomal.loc <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene") [genes]
promoter.gene <- getPromoterSeq(chromosomal.loc, Mmusculus, upstream=2000, downstream=1000)

###ANdrogen weight matrices
weights.ADR3 <- rbind(
c(1,0,1,1,0,1,0.25,0.25,0.25,1,0,1,1,0,1), #A
c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,0), #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0), #G
c(0,0,0,0,0,0,0.25,0.25,0.25,0,0,0,0,0,0)) #T
rownames(weights.ADR3) <- c("A","C","G","T")
weights.ARE2 <- rbind(c(1,0,0,1,0,0,0.25,0.25,0.25,1,0,1,1,0,0), #A
c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,1), #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0), #G
c(0,0,1,0,0,1,0.25,0.25,0.25,0,0,0,0,0,0)) #T
rownames(weights.ARE2) <- c("A","C","G","T")
weights.ARE <- rbind(
c(1,0,1,1,0,1,0.25,0.25,0.25,0,0,0,1,0,0), #A
c(0,0,0,0,1,0,0.25,0.25,0.25,1,0,0,0,1,0), #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0), #G
c(0,0,0,0,0,0,0.25,0.25,0.25,0,0,1,0,0,1)) #T
rownames(weights.ARE) <- c("A","C","G","T")

weights.IR3 <- rbind(
c(1,0,1,1,0,1,0.25,0.25,0.25,0,0,0,0,0,0), #A
c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,0), #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,1), #G
c(0,0,0,0,0,0,0.25,0.25,0.25,1,0,1,1,0,0)) #T
rownames(weights.IR3) <- c("A","C","G","T")


###
# "simple" estrogen T00259 ER-alpha; Species: mouse, Mus musculus.http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=m00191
#weights.ER <- rbind(
#c(0.25,0.25,0.6,0.3,0.15,0.25,0.15,0.55,0.25,0.25,0.25,0,0.05,1,0,0.05,0.1,0.25,0.25),	#A
#c(0.25,0.25,0.2,0.05,0.1,0.25,0.5,0.1,0.25,0.25,0.25,0.05,0,0,0.9,0.9,0.35,0.25,0.25), 	#C
#c(0.25,0.25,0.15,0.55,0.55,0.25,0.15,0.2,0.25,0.25,0.25,0,0.95,0,0.1,0.05,0.1,0.25,0.25), 	#G
#c(0.25,0.25,0.05,0.1,0.2,0.25,0.2,0.15,0.25,0.25,0.25,0.95,0,0,0,0,0.45,0.25,0.25))	#T
#rownames(weights.ER) <- c("A","C","G","T")
# "complex" ER, diff are the N %
weights.ER <- rbind(
c(0.3,0.2,0.6,0.3,0.15,0.15,0.15,0.55,0.2,0.15,0.1,0,0.05,1,0,0.05,0.1,0.35,0.4),	#A
c(0.2,0.45,0.2,0.05,0.1,0.15,0.5,0.1,0.45,0.3,0.15,0.05,0,0,0.9,0.9,0.35,0.05,0.25),	#C
c(0.2,0.15,0.15,0.55,0.55,0.2,0.15,0.2,0.15,0.15,0.45,0,0.95,0,0.1,0.05,0.1,0.4,0.15),	#G
c(0.3,0.2,0.05,0.1,0.2,0.5,0.2,0.15,0.2,0.4,0.3,0.95,0,0,0,0,0.45,0.2,0.2))	#T
rownames(weights.ER) <- c("A","C","G","T")

# "simple" Glucocorticoid receptor, http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=M00192
#weights.GR <- rbind(
#c(0.25,0.25,0.25,0.25,0.25,0.25,0.13,0.25,0.25,0.26,0.25,0,0,0.11,0.25,0.05,0.24,0.25,0.25),	#A
#c(0.25,0.25,0.25,0.25,0.25,0.25,0.58,0.25,0.25,0.08,0.25,0,0,0.03,0.25,0.79,0.13,0.25,0.25), 	#C
#c(0.25,0.25,0.25,0.25,0.25,0.25,0.16,0.25,0.25,0.08,0.25,0,1,0,0.25,0.13,0,0.25,0.25), 	#G
#c(0.25,0.25,0.25,0.25,0.25,0.25,0.13,0.25,0.25,0.58,0.25,1,0,0.86,0.25,0.03,0.63,0.25,0.25))	#T
#rownames(weights.GR) <- c("A","C","G","T")
# "complex" GR, diff are the N %
weights.GR <- rbind(
c(0.24,0.21,0.21,0.24,0.24,0.37,0.13,0.39,0.32,0.26,0.21,0,0,0.11,0.10,0.05,0.24,0.32,0.20),	#A
c(0.15,0.21,0.11,0.21,0.24,0.24,0.58,0.21,0.34,0.08,0.26,0,0,0.02,0.24,0.79,0.13,0.21,0.24),	#C
c(0.24,0.37,0.47,0.45,0.18,0.10,0.16,0.16,0.10,0.08,0.39,0,1,0,0.16,0.13,0,0.15,0.32),	#G
c(0.37,0.21,0.21,0.10,0.34,0.29,0.13,0.24,0.24,0.58,0.14,1,0,0.87,0.50,0.03,0.63,0.32,0.24))	#T
rownames(weights.GR) <- c("A","C","G","T")

# ERE estr resp elem. (Chin-Yo Lin cartography)
weights.ERE <- rbind(
c(0,0,0,0,1,0.25,0.25,0.25,0,0,1,0,0),	#A
c(0,0,0,1,0,0.25,0.25,0.25,0,0,0,1,1),	#C
c(1,1,0,0,0,0.25,0.25,0.25,0,1,0,0,0),	#G
c(0,0,1,0,0,0.25,0.25,0.25,1,0,0,0,0))	#T
rownames(weights.ERE) <- c("A","C","G","T")

# ESR1 MA0112.1 JASPAR
weights.ESR1 <- rbind(
c(0.11,0.11,0.78,0.22,0,0,0,0.67,0.11,0.22,0.33,0.11,0.11,0.56,0,0,0.22,0.33),	#A
c(0.56,0.56,0.11,0,0,0,0.78,0,0.78,0.56,0.22,0.11,0,0.11,0.89,1,0.45,0.45),	#C
c(0.11,0.11,0.11,0.78,1,0,0.22,0.22,0.11,0.11,0.45,0.11,0.89,0.33,0,0,0,0.11),	#G
c(0.22,0.22,0,0,0,1,0,0.11,0,0.11,0,0.67,0,0,0.11,0,0.33,0.11))	#T
rownames(weights.ESR1) <- c("A","C","G","T")

wmatrices <- list(weights.ADR3, weights.ARE2, weights.ARE, weights.IR3, weights.ER, weights.GR, weights.ERE,weights.ESR1)
names(wmatrices) <- c("ADR3","ARE2","ARE","IR3", "ER", "GR", "ERE","ESR1")
output <- NULL
for(x in 1:length(wmatrices)){
seqLogo(wmatrices[[x]])
TFBSjaspar <- round(100 *wmatrices[[x]])
allhits <- NULL
for(gene in names(promoter.gene)){
hits <- matchPWM(TFBSjaspar, unlist(promoter.gene[gene])[[1]], "90%")
allhits <- c(allhits, length(hits))
}
output <- cbind(output, allhits)
}
rownames(output) <- names(promoter.gene)
colnames(output) <- names(wmatrices)
#library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") # Biomart for mouse genes
fullnames <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), filters = c("ensembl_gene_id"), values = rownames(output), mart = bio.mart)
output <- output[match(fullnames[,"ensembl_gene_id"], rownames(output)),]
write.table(cbind(fullnames[,"mgi_symbol"], output), file="TFBSmatches_PaperGenes_AREundcomplexERundcomplexGRundEREundESR1.90prz.txt",sep="\t",quote=FALSE)


### paper genes 
genes <- c("ENSMUSG00000020102","ENSMUSG00000022853","ENSMUSG00000020777","ENSMUSG00000024286","ENSMUSG00000024066","ENSMUSG00000009646","ENSMUSG00000067369","ENSMUSG00000031446","ENSMUSG00000028127","ENSMUSG00000031996","ENSMUSG00000038982","ENSMUSG00000026072","ENSMUSG00000049044",
"ENSMUSG00000057236","ENSMUSG00000009281","ENSMUSG00000023031","ENSMUSG00000002014","ENSMUSG00000028405","ENSMUSG00000045636","ENSMUSG00000000751","ENSMUSG00000057841","ENSMUSG00000038495","ENSMUSG00000078817","ENSMUSG00000031320","ENSMUSG00000038521","ENSMUSG00000029370",
"ENSMUSG00000002996","ENSMUSG00000033174","ENSMUSG00000065979","ENSMUSG00000025981","ENSMUSG00000081043","ENSMUSG00000083563","ENSMUSG00000010025","ENSMUSG00000027332","ENSMUSG00000032067","ENSMUSG00000009907","ENSMUSG00000028158","ENSMUSG00000022707","ENSMUSG00000026827",
"ENSMUSG00000026638","ENSMUSG00000051147","ENSMUSG00000035093","ENSMUSG00000050069","ENSMUSG00000037287","ENSMUSG00000031024","ENSMUSG00000025289","ENSMUSG00000032310","ENSMUSG00000020553","ENSMUSG00000025153","ENSMUSG00000038298","ENSMUSG00000029998","ENSMUSG00000055172",
"ENSMUSG00000064341","ENSMUSG00000015112","ENSMUSG00000027984","ENSMUSG00000020986","ENSMUSG00000046598","ENSMUSG00000025745","ENSMUSG00000001783","ENSMUSG00000047676","ENSMUSG00000078517","ENSMUSG00000021559","ENSMUSG00000028024","ENSMUSG00000090150","ENSMUSG00000068877",
"ENSMUSG00000021236","ENSMUSG00000072849","ENSMUSG00000025453","ENSMUSG00000021136","ENSMUSG00000032418","ENSMUSG00000023043","ENSMUSG00000027439","ENSMUSG00000028980","ENSMUSG00000033107","ENSMUSG00000028527","ENSMUSG00000058638","ENSMUSG00000045983","ENSMUSG00000016664",
"ENSMUSG00000026853","ENSMUSG00000019935","ENSMUSG00000020496","ENSMUSG00000029767","ENSMUSG00000037053","ENSMUSG00000038155","ENSMUSG00000074272","ENSMUSG00000028760","ENSMUSG00000022790","ENSMUSG00000058301","ENSMUSG00000016024","ENSMUSG00000032883","ENSMUSG00000032193",
"ENSMUSG00000051065","ENSMUSG00000041220","ENSMUSG00000029672","ENSMUSG00000071646","ENSMUSG00000027309","ENSMUSG00000009905","ENSMUSG00000032802","ENSMUSG00000020576","ENSMUSG00000063354","ENSMUSG00000022893","ENSMUSG00000055782","ENSMUSG00000061313","ENSMUSG00000038591",
"ENSMUSG00000021670","ENSMUSG00000041957","ENSMUSG00000032207","ENSMUSG00000020534","ENSMUSG00000019122","ENSMUSG00000070939","ENSMUSG00000078672","ENSMUSG00000024887","ENSMUSG00000003970","ENSMUSG00000038764","ENSMUSG00000035948","ENSMUSG00000054263","ENSMUSG00000030934",
"ENSMUSG00000022512","ENSMUSG00000030435","ENSMUSG00000084309","ENSMUSG00000029556","ENSMUSG00000020122","ENSMUSG00000020521","ENSMUSG00000025885","ENSMUSG00000010651","ENSMUSG00000067279","ENSMUSG00000022821","ENSMUSG00000034371","ENSMUSG00000084983","ENSMUSG00000022246",
"ENSMUSG00000024986","ENSMUSG00000096188","ENSMUSG00000019951","ENSMUSG00000028223","ENSMUSG00000000303","ENSMUSG00000060803","ENSMUSG00000020300","ENSMUSG00000061292","ENSMUSG00000032702","ENSMUSG00000038068","ENSMUSG00000014542","ENSMUSG00000000826","ENSMUSG00000090136",
"ENSMUSG00000058558","ENSMUSG00000055116","ENSMUSG00000097657","ENSMUSG00000040713","ENSMUSG00000035686","ENSMUSG00000064225","ENSMUSG00000063683","ENSMUSG00000037573","ENSMUSG00000027312","ENSMUSG00000024507","ENSMUSG00000032232","ENSMUSG00000035473","ENSMUSG00000047514",
"ENSMUSG00000039115","ENSMUSG00000090555","ENSMUSG00000029153","ENSMUSG00000032194","ENSMUSG00000033860","ENSMUSG00000027472","ENSMUSG00000022032","ENSMUSG00000063882","ENSMUSG00000045180","ENSMUSG00000050730","ENSMUSG00000047446","ENSMUSG00000032235","ENSMUSG00000082644",
"ENSMUSG00000089991","ENSMUSG00000067144","ENSMUSG00000035674","ENSMUSG00000047407","ENSMUSG00000062181","ENSMUSG00000025036","ENSMUSG00000026675","ENSMUSG00000027642","ENSMUSG00000045092","ENSMUSG00000095042","ENSMUSG00000018425","ENSMUSG00000037058","ENSMUSG00000020390",
"ENSMUSG00000049164","ENSMUSG00000035133","ENSMUSG00000081049","ENSMUSG00000004552","ENSMUSG00000035203","ENSMUSG00000046959","ENSMUSG00000021411","ENSMUSG00000020865","ENSMUSG00000027762","ENSMUSG00000034211","ENSMUSG00000029076","ENSMUSG00000037712","ENSMUSG00000025810",
"ENSMUSG00000006273","ENSMUSG00000019795","ENSMUSG00000035875","ENSMUSG00000052392","ENSMUSG00000017009","ENSMUSG00000046352","ENSMUSG00000009378","ENSMUSG00000026874","ENSMUSG00000087370","ENSMUSG00000058291","ENSMUSG00000015305","ENSMUSG00000032808","ENSMUSG00000031604",
"ENSMUSG00000025809","ENSMUSG00000034667","ENSMUSG00000024978","ENSMUSG00000021102","ENSMUSG00000030744","ENSMUSG00000060429","ENSMUSG00000026848","ENSMUSG00000058838","ENSMUSG00000001666","ENSMUSG00000025950","ENSMUSG00000039952","ENSMUSG00000071655","ENSMUSG00000053536",
"ENSMUSG00000032116","ENSMUSG00000096449","ENSMUSG00000021226","ENSMUSG00000055093","ENSMUSG00000037818","ENSMUSG00000020672","ENSMUSG00000058624","ENSMUSG00000037905","ENSMUSG00000041293","ENSMUSG00000080738","ENSMUSG00000005951","ENSMUSG00000021497","ENSMUSG00000052656",
"ENSMUSG00000055730","ENSMUSG00000043716","ENSMUSG00000041236","ENSMUSG00000034528","ENSMUSG00000015357","ENSMUSG00000003429","ENSMUSG00000027397","ENSMUSG00000055737","ENSMUSG00000041567","ENSMUSG00000028653","ENSMUSG00000034926","ENSMUSG00000047492","ENSMUSG00000047675",
"ENSMUSG00000033308","ENSMUSG00000020257","ENSMUSG00000035129","ENSMUSG00000029656","ENSMUSG00000016831","ENSMUSG00000084038","ENSMUSG00000068086","ENSMUSG00000030513","ENSMUSG00000041390","ENSMUSG00000029686","ENSMUSG00000049106","ENSMUSG00000028957","ENSMUSG00000029173",
"ENSMUSG00000006445","ENSMUSG00000002332","ENSMUSG00000062382","ENSMUSG00000009614","ENSMUSG00000005225","ENSMUSG00000041528","ENSMUSG00000040466","ENSMUSG00000033720","ENSMUSG00000025702","ENSMUSG00000072949","ENSMUSG00000036905","ENSMUSG00000031451","ENSMUSG00000001604",
"ENSMUSG00000031371","ENSMUSG00000060317","ENSMUSG00000034780","ENSMUSG00000063558","ENSMUSG00000033382","ENSMUSG00000034837","ENSMUSG00000049502","ENSMUSG00000022280","ENSMUSG00000026435","ENSMUSG00000005615","ENSMUSG00000044345","ENSMUSG00000028199","ENSMUSG00000049323",
"ENSMUSG00000024132","ENSMUSG00000024335","ENSMUSG00000020572","ENSMUSG00000021751","ENSMUSG00000056131","ENSMUSG00000010601","ENSMUSG00000006736","ENSMUSG00000024608","ENSMUSG00000041130","ENSMUSG00000032399","ENSMUSG00000053931","ENSMUSG00000049577","ENSMUSG00000050079",
"ENSMUSG00000025003","ENSMUSG00000029787","ENSMUSG00000039428","ENSMUSG00000010517","ENSMUSG00000091561","ENSMUSG00000038145","ENSMUSG00000030678","ENSMUSG00000022629","ENSMUSG00000022149","ENSMUSG00000039001","ENSMUSG00000024391","ENSMUSG00000000605","ENSMUSG00000066154",
"ENSMUSG00000054452","ENSMUSG00000033916","ENSMUSG00000054072","ENSMUSG00000010663","ENSMUSG00000020590","ENSMUSG00000027710","ENSMUSG00000080832","ENSMUSG00000033629","ENSMUSG00000040029","ENSMUSG00000032536","ENSMUSG00000078920","ENSMUSG00000003585","ENSMUSG00000038206",
"ENSMUSG00000037075","ENSMUSG00000002881","ENSMUSG00000022766","ENSMUSG00000053460","ENSMUSG00000033400","ENSMUSG00000038092","ENSMUSG00000015522","ENSMUSG00000020064","ENSMUSG00000028826","ENSMUSG00000055866","ENSMUSG00000063684","ENSMUSG00000027809","ENSMUSG00000054823",
"ENSMUSG00000067274","ENSMUSG00000035776","ENSMUSG00000047905","ENSMUSG00000011831","ENSMUSG00000028063","ENSMUSG00000070594","ENSMUSG00000027774","ENSMUSG00000084319","ENSMUSG00000005936","ENSMUSG00000032238","ENSMUSG00000085867","ENSMUSG00000025936","ENSMUSG00000098470",
"ENSMUSG00000026618","ENSMUSG00000071650","ENSMUSG00000026112","ENSMUSG00000017929","ENSMUSG00000080776","ENSMUSG00000030556","ENSMUSG00000010122","ENSMUSG00000028961","ENSMUSG00000040181","ENSMUSG00000064373","ENSMUSG00000028284","ENSMUSG00000004040","ENSMUSG00000066072",
"ENSMUSG00000017132","ENSMUSG00000022389","ENSMUSG00000061778","ENSMUSG00000041660","ENSMUSG00000016382","ENSMUSG00000072082","ENSMUSG00000027678","ENSMUSG00000032199","ENSMUSG00000026811","ENSMUSG00000003477","ENSMUSG00000061838","ENSMUSG00000024283","ENSMUSG00000054648",
"ENSMUSG00000090306","ENSMUSG00000019188","ENSMUSG00000060938","ENSMUSG00000057637","ENSMUSG00000018846","ENSMUSG00000021823","ENSMUSG00000028792","ENSMUSG00000034761","ENSMUSG00000026623","ENSMUSG00000024621","ENSMUSG00000023030","ENSMUSG00000063760","ENSMUSG00000025854",
"ENSMUSG00000022956","ENSMUSG00000032198","ENSMUSG00000029238","ENSMUSG00000021973","ENSMUSG00000071204","ENSMUSG00000031605","ENSMUSG00000030237","ENSMUSG00000024187","ENSMUSG00000029763","ENSMUSG00000036368","ENSMUSG00000040017","ENSMUSG00000035086","ENSMUSG00000007656",
"ENSMUSG00000026991","ENSMUSG00000017760","ENSMUSG00000091549","ENSMUSG00000046691","ENSMUSG00000069682","ENSMUSG00000033377","ENSMUSG00000033955","ENSMUSG00000054455","ENSMUSG00000020115","ENSMUSG00000020021","ENSMUSG00000040127","ENSMUSG00000021871","ENSMUSG00000067235",
"ENSMUSG00000031349","ENSMUSG00000008668","ENSMUSG00000050708","ENSMUSG00000020166","ENSMUSG00000029063","ENSMUSG00000039046","ENSMUSG00000020694","ENSMUSG00000029833","ENSMUSG00000049555","ENSMUSG00000021327","ENSMUSG00000041895","ENSMUSG00000046532","ENSMUSG00000036138",
"ENSMUSG00000028051","ENSMUSG00000059291","ENSMUSG00000077637","ENSMUSG00000039519","ENSMUSG00000045868","ENSMUSG00000020429","ENSMUSG00000039349","ENSMUSG00000026295","ENSMUSG00000025260","ENSMUSG00000073640","ENSMUSG00000049775","ENSMUSG00000023186","ENSMUSG00000016252",
"ENSMUSG00000020017","ENSMUSG00000039967","ENSMUSG00000003228","ENSMUSG00000022957","ENSMUSG00000079470","ENSMUSG00000032409","ENSMUSG00000081824","ENSMUSG00000021585","ENSMUSG00000025584","ENSMUSG00000074882","ENSMUSG00000025937","ENSMUSG00000026701","ENSMUSG00000024781",
"ENSMUSG00000022967","ENSMUSG00000048222","ENSMUSG00000021273","ENSMUSG00000048600","ENSMUSG00000033326","ENSMUSG00000064437","ENSMUSG00000024843","ENSMUSG00000032602","ENSMUSG00000031672","ENSMUSG00000026986","ENSMUSG00000015189","ENSMUSG00000039000","ENSMUSG00000071551",
"ENSMUSG00000053897","ENSMUSG00000031844","ENSMUSG00000033688","ENSMUSG00000037953","ENSMUSG00000033965","ENSMUSG00000039275","ENSMUSG00000033285","ENSMUSG00000038612","ENSMUSG00000020460","ENSMUSG00000038696","ENSMUSG00000024603","ENSMUSG00000020653","ENSMUSG00000024566",
"ENSMUSG00000024493","ENSMUSG00000004268","ENSMUSG00000066415","ENSMUSG00000022181","ENSMUSG00000010097","ENSMUSG00000020034","ENSMUSG00000096674","ENSMUSG00000026272","ENSMUSG00000070565","ENSMUSG00000027075","ENSMUSG00000079362","ENSMUSG00000057666","ENSMUSG00000031818",
"ENSMUSG00000047963","ENSMUSG00000053617","ENSMUSG00000006529","ENSMUSG00000060126","ENSMUSG00000046311","ENSMUSG00000046324","ENSMUSG00000002985","ENSMUSG00000034781","ENSMUSG00000026687","ENSMUSG00000073460","ENSMUSG00000030108","ENSMUSG00000056579","ENSMUSG00000020265",
"ENSMUSG00000024982","ENSMUSG00000025014","ENSMUSG00000031645","ENSMUSG00000016028","ENSMUSG00000079509","ENSMUSG00000052539","ENSMUSG00000037025","ENSMUSG00000024236","ENSMUSG00000029246","ENSMUSG00000069972","ENSMUSG00000039159","ENSMUSG00000094248","ENSMUSG00000036611",
"ENSMUSG00000020917","ENSMUSG00000052974","ENSMUSG00000025533","ENSMUSG00000023019","ENSMUSG00000038642","ENSMUSG00000030869","ENSMUSG00000020538","ENSMUSG00000027253","ENSMUSG00000031378","ENSMUSG00000028081","ENSMUSG00000021215","ENSMUSG00000057342","ENSMUSG00000045294",
"ENSMUSG00000048039","ENSMUSG00000030064","ENSMUSG00000036151","ENSMUSG00000020123","ENSMUSG00000036896","ENSMUSG00000078965","ENSMUSG00000052914","ENSMUSG00000056268","ENSMUSG00000056973","ENSMUSG00000038876","ENSMUSG00000069805","ENSMUSG00000039047","ENSMUSG00000029020",
"ENSMUSG00000021999","ENSMUSG00000021366","ENSMUSG00000031271","ENSMUSG00000016541","ENSMUSG00000035399","ENSMUSG00000028211","ENSMUSG00000029482","ENSMUSG00000005514","ENSMUSG00000009927","ENSMUSG00000026858","ENSMUSG00000003518","ENSMUSG00000031441","ENSMUSG00000032724",
"ENSMUSG00000025317","ENSMUSG00000097245","ENSMUSG00000078606","ENSMUSG00000038776","ENSMUSG00000030946","ENSMUSG00000039197","ENSMUSG00000021993","ENSMUSG00000066026","ENSMUSG00000042228","ENSMUSG00000067547","ENSMUSG00000024863","ENSMUSG00000042363","ENSMUSG00000059810",
"ENSMUSG00000034993","ENSMUSG00000030787","ENSMUSG00000032383","ENSMUSG00000026368","ENSMUSG00000006494","ENSMUSG00000023057","ENSMUSG00000020571","ENSMUSG00000024665","ENSMUSG00000039512","ENSMUSG00000037470","ENSMUSG00000024424","ENSMUSG00000032018","ENSMUSG00000004730",
"ENSMUSG00000018547","ENSMUSG00000028715","ENSMUSG00000050777","ENSMUSG00000057367","ENSMUSG00000001630","ENSMUSG00000027951","ENSMUSG00000031561","ENSMUSG00000065852","ENSMUSG00000025059","ENSMUSG00000089407","ENSMUSG00000005968","ENSMUSG00000029630","ENSMUSG00000022848",
"ENSMUSG00000039221","ENSMUSG00000038768","ENSMUSG00000024661","ENSMUSG00000059481","ENSMUSG00000051483","ENSMUSG00000029385","ENSMUSG00000081221","ENSMUSG00000040564","ENSMUSG00000019718","ENSMUSG00000082100","ENSMUSG00000070730","ENSMUSG00000059908","ENSMUSG00000042364",
"ENSMUSG00000000056","ENSMUSG00000001100","ENSMUSG00000033634","ENSMUSG00000024052","ENSMUSG00000078636","ENSMUSG00000029775","ENSMUSG00000015971","ENSMUSG00000020532","ENSMUSG00000012405","ENSMUSG00000022663","ENSMUSG00000032735","ENSMUSG00000045767","ENSMUSG00000027875",
"ENSMUSG00000058230","ENSMUSG00000045098","ENSMUSG00000079343","ENSMUSG00000021033","ENSMUSG00000078193","ENSMUSG00000045896","ENSMUSG00000021917","ENSMUSG00000027423","ENSMUSG00000036887","ENSMUSG00000021583","ENSMUSG00000020375","ENSMUSG00000068739","ENSMUSG00000002147",
"ENSMUSG00000022272","ENSMUSG00000010936","ENSMUSG00000052040","ENSMUSG00000027765","ENSMUSG00000031444","ENSMUSG00000038224","ENSMUSG00000026849","ENSMUSG00000039361","ENSMUSG00000032078","ENSMUSG00000020910","ENSMUSG00000029247","ENSMUSG00000031166","ENSMUSG00000002265",
"ENSMUSG00000031157","ENSMUSG00000071644","ENSMUSG00000022844","ENSMUSG00000020700","ENSMUSG00000016534","ENSMUSG00000019189","ENSMUSG00000078300","ENSMUSG00000024772","ENSMUSG00000056537","ENSMUSG00000039682","ENSMUSG00000019302","ENSMUSG00000020988","ENSMUSG00000029198",
"ENSMUSG00000031443","ENSMUSG00000060807","ENSMUSG00000036503","ENSMUSG00000095937","ENSMUSG00000067870","ENSMUSG00000031644","ENSMUSG00000057228","ENSMUSG00000035769","ENSMUSG00000051375","ENSMUSG00000037876","ENSMUSG00000025935","ENSMUSG00000021124","ENSMUSG00000026473",
"ENSMUSG00000021814","ENSMUSG00000015755","ENSMUSG00000009621","ENSMUSG00000020009","ENSMUSG00000022514","ENSMUSG00000075268","ENSMUSG00000070644","ENSMUSG00000030545","ENSMUSG00000005575","ENSMUSG00000047230","ENSMUSG00000019851","ENSMUSG00000075706","ENSMUSG00000016940",
"ENSMUSG00000025813","ENSMUSG00000081232","ENSMUSG00000024610","ENSMUSG00000019877","ENSMUSG00000078087","ENSMUSG00000032548","ENSMUSG00000001665","ENSMUSG00000094708","ENSMUSG00000010048","ENSMUSG00000033416","ENSMUSG00000026568","ENSMUSG00000040888","ENSMUSG00000088609",
"ENSMUSG00000015970","ENSMUSG00000020620","ENSMUSG00000020530","ENSMUSG00000028646","ENSMUSG00000068134","ENSMUSG00000082051","ENSMUSG00000020087","ENSMUSG00000039529","ENSMUSG00000021210","ENSMUSG00000062825","ENSMUSG00000060450","ENSMUSG00000035027","ENSMUSG00000051329",
"ENSMUSG00000074582","ENSMUSG00000053846","ENSMUSG00000039114","ENSMUSG00000036918","ENSMUSG00000003623","ENSMUSG00000038806","ENSMUSG00000027546","ENSMUSG00000033943","ENSMUSG00000021576","ENSMUSG00000035031","ENSMUSG00000024085","ENSMUSG00000032407","ENSMUSG00000035936",
"ENSMUSG00000046330","ENSMUSG00000021417","ENSMUSG00000041237","ENSMUSG00000021281","ENSMUSG00000054422","ENSMUSG00000043131","ENSMUSG00000034714","ENSMUSG00000093930","ENSMUSG00000006333","ENSMUSG00000036104","ENSMUSG00000060301","ENSMUSG00000063590","ENSMUSG00000082062",
"ENSMUSG00000022771","ENSMUSG00000032420","ENSMUSG00000057103","ENSMUSG00000027601","ENSMUSG00000030527","ENSMUSG00000033488","ENSMUSG00000004897","ENSMUSG00000036112","ENSMUSG00000032249","ENSMUSG00000095478","ENSMUSG00000025792","ENSMUSG00000071076","ENSMUSG00000021690",
"ENSMUSG00000045055","ENSMUSG00000006134","ENSMUSG00000058997","ENSMUSG00000078686","ENSMUSG00000031885","ENSMUSG00000055991","ENSMUSG00000025396","ENSMUSG00000067071","ENSMUSG00000001323","ENSMUSG00000056666","ENSMUSG00000039157","ENSMUSG00000021096","ENSMUSG00000066880",
"ENSMUSG00000024845","ENSMUSG00000074240","ENSMUSG00000037798","ENSMUSG00000022946","ENSMUSG00000067212","ENSMUSG00000023861","ENSMUSG00000097679","ENSMUSG00000022037","ENSMUSG00000044350","ENSMUSG00000006673","ENSMUSG00000052155","ENSMUSG00000070394","ENSMUSG00000073842",
"ENSMUSG00000069456","ENSMUSG00000047246","ENSMUSG00000029864","ENSMUSG00000020463","ENSMUSG00000025162","ENSMUSG00000029614","ENSMUSG00000091264","ENSMUSG00000026348","ENSMUSG00000033707","ENSMUSG00000027801","ENSMUSG00000020010","ENSMUSG00000023079","ENSMUSG00000021361",
"ENSMUSG00000006717","ENSMUSG00000022791","ENSMUSG00000003053","ENSMUSG00000049804","ENSMUSG00000083773","ENSMUSG00000038459","ENSMUSG00000015247","ENSMUSG00000027357","ENSMUSG00000039844","ENSMUSG00000037286","ENSMUSG00000022684","ENSMUSG00000039886","ENSMUSG00000061175",
"ENSMUSG00000053870","ENSMUSG00000046364","ENSMUSG00000027690","ENSMUSG00000028381","ENSMUSG00000037519","ENSMUSG00000027890","ENSMUSG00000071267","ENSMUSG00000026889","ENSMUSG00000028076","ENSMUSG00000019818","ENSMUSG00000015890","ENSMUSG00000030016","ENSMUSG00000044068",
"ENSMUSG00000056204","ENSMUSG00000044477","ENSMUSG00000038267","ENSMUSG00000031173","ENSMUSG00000038429","ENSMUSG00000027533","ENSMUSG00000022210","ENSMUSG00000051391","ENSMUSG00000087141","ENSMUSG00000095361","ENSMUSG00000018427","ENSMUSG00000002265","ENSMUSG00000004562",
"ENSMUSG00000006289","ENSMUSG00000016024","ENSMUSG00000019188","ENSMUSG00000022789","ENSMUSG00000025066","ENSMUSG00000025439","ENSMUSG00000026659","ENSMUSG00000028715","ENSMUSG00000029314","ENSMUSG00000030131","ENSMUSG00000032271","ENSMUSG00000037440","ENSMUSG00000041278",
"ENSMUSG00000041567","ENSMUSG00000042248","ENSMUSG00000044068","ENSMUSG00000046667","ENSMUSG00000048758","ENSMUSG00000057037","ENSMUSG00000057074","ENSMUSG00000058126","ENSMUSG00000061232","ENSMUSG00000062054","ENSMUSG00000063232","ENSMUSG00000063590","ENSMUSG00000064193",
"ENSMUSG00000066366","ENSMUSG00000069274","ENSMUSG00000071177","ENSMUSG00000073418","ENSMUSG00000074272","ENSMUSG00000074882","ENSMUSG00000079427","ENSMUSG00000081207","ENSMUSG00000090555","ENSMUSG00000091537","ENSMUSG00000092021","ENSMUSG00000093753","ENSMUSG00000054676",
"ENSMUSG00000097354","ENSMUSG00000050541","ENSMUSG00000029314","ENSMUSG00000026542","ENSMUSG00000071176","ENSMUSG00000004562","ENSMUSG00000028008","ENSMUSG00000022533","ENSMUSG00000028457","ENSMUSG00000073418","ENSMUSG00000022181","ENSMUSG00000102095","ENSMUSG00000032803",
"ENSMUSG00000074272","ENSMUSG00000037443","ENSMUSG00000057074","ENSMUSG00000057037","ENSMUSG00000025439","ENSMUSG00000016756","ENSMUSG00000030560","ENSMUSG00000042248","ENSMUSG00000074882","ENSMUSG00000021259","ENSMUSG00000028715","ENSMUSG00000044748","ENSMUSG00000022789",
"ENSMUSG00000026659","ENSMUSG00000048376","ENSMUSG00000072568","ENSMUSG00000028218","ENSMUSG00000092021","ENSMUSG00000081207","ENSMUSG00000079427","ENSMUSG00000064193","ENSMUSG00000090555","ENSMUSG00000025190","ENSMUSG00000019188","ENSMUSG00000061232","ENSMUSG00000030532",
"ENSMUSG00000042770","ENSMUSG00000024986","ENSMUSG00000069274","ENSMUSG00000037234","ENSMUSG00000062054","ENSMUSG00000045294","ENSMUSG00000035578","ENSMUSG00000016024","ENSMUSG00000020593","ENSMUSG00000027286","ENSMUSG00000069516","ENSMUSG00000054387","ENSMUSG00000009376",
"ENSMUSG00000008540","ENSMUSG00000027820","ENSMUSG00000039395","ENSMUSG00000034729","ENSMUSG00000030131","ENSMUSG00000037847","ENSMUSG00000032271","ENSMUSG00000028851","ENSMUSG00000006289","ENSMUSG00000037348","ENSMUSG00000056851","ENSMUSG00000031592","ENSMUSG00000020553",
"ENSMUSG00000002265","ENSMUSG00000005161","ENSMUSG00000025314","ENSMUSG00000046667","ENSMUSG00000040782","ENSMUSG00000048758","ENSMUSG00000063253","ENSMUSG00000022772","ENSMUSG00000063232","ENSMUSG00000041567","ENSMUSG00000066366","ENSMUSG00000071177","ENSMUSG00000058207",
"ENSMUSG00000025066","ENSMUSG00000036083","ENSMUSG00000033147","ENSMUSG00000063590","ENSMUSG00000032602","ENSMUSG00000027227","ENSMUSG00000091537","ENSMUSG00000058126","ENSMUSG00000041278","ENSMUSG00000031577","ENSMUSG00000023286","ENSMUSG00000041231","ENSMUSG00000059534",
"ENSMUSG00000028452","ENSMUSG00000037440","ENSMUSG00000071281","ENSMUSG00000044068")


RPKM[which(RPKM$ensembl_gene_id %in% genes),1]
length(genes)

### estrogen ER http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=m00191


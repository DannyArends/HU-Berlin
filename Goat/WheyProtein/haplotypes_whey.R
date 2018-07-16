library("haplo.stats")

group <- as.factor(c(rep("Nilotic", 5), rep("Desert", 5),rep("Taggar", 5),rep("Nubian", 5)))
setwd("D:/Edrive/Goat/DNA/WheyHaplotypes")

lalba_upstream <- read.table("lalba_upstream.txt",sep="\t", colClasses="character")
lalba_exonic <- read.table("lalba_exonic.txt",sep="\t", colClasses="character")
lgb_upstream <- read.table("lgb_upstream.txt",sep="\t", colClasses="character")
lgb_exonic <- read.table("lgb_exonic.txt",sep="\t", colClasses="character")

labels_lalba_upstream = c("rs642745519", "rs657752748", "rs672221003", "rs660607813", "rs646115281", "rs671985793", "rs649340231")
labels_lalba_exonic = c("rs660007102", "rs650426572", "5:30813099A>G", "5:30813117A>G", "5:30813123A>G", "rs663851286", "rs657915405", "rs641559728", "rs664225585")

labels_lgb_upstream = c("11:102720926A>G", "rs646608188", "rs659709384", "rs659299918", "rs635615192", "11:102722042G>T", "11:102722111C>T", "11:102722160G>A")
labels_lgb_exonic = c("11:102723010G>A", "rs655422073", "rs654002978", "rs651951335", "rs666423193")

res_lalba_upstream <- haplo.group(group, lalba_upstream, labels_lalba_upstream, haplo.em.control(max.iter = 100, min.posterior = 1e-06, verbose=TRUE))
res_lalba_exonic <- haplo.group(group, lalba_exonic, labels_lalba_exonic, haplo.em.control(max.iter = 100, min.posterior = 1e-06, verbose=TRUE))

res_lgb_upstream <- haplo.group(group, lgb_upstream, labels_lgb_upstream, haplo.em.control(max.iter = 100, min.posterior = 1e-06, verbose=TRUE))
res_lgb_exonic <- haplo.group(group, lgb_exonic, labels_lgb_exonic, haplo.em.control(max.iter = 100, min.posterior = 1e-06, verbose=TRUE))

#lalba_upstream
lalba_upstream_haplo = res_lalba_upstream[[1]]
colnames(lalba_upstream_haplo) <- gsub("group=", "", colnames(lalba_upstream_haplo))
lalba_upstream_haplo[,"Total"] =  round(lalba_upstream_haplo[,"Total"], 3)
lalba_upstream_haplo[,"Desert"] =  round(lalba_upstream_haplo[,"Desert"], 3)
lalba_upstream_haplo[,"Taggar"] =  round(lalba_upstream_haplo[,"Taggar"], 3)
lalba_upstream_haplo[,"Nubian"] =  round(lalba_upstream_haplo[,"Nubian"], 3)
lalba_upstream_haplo[,"Nilotic"] =  round(lalba_upstream_haplo[,"Nilotic"], 3)

notPresent = which(apply(lalba_upstream_haplo[,c("Total", as.character(unique(group)))],1,sum,na.rm=T) == 0)
lalba_upstream_haplo = lalba_upstream_haplo[-notPresent,]
lalba_upstream_haplo = lalba_upstream_haplo[order(lalba_upstream_haplo[, "Total"], decreasing = TRUE),]
rownames(lalba_upstream_haplo) <- paste("Haplotype", LETTERS[1:nrow(lalba_upstream_haplo)])

write.table(lalba_upstream_haplo, "lalba_haplotypes_upstream.txt", sep = "\t", quote=FALSE, na="")

#lalba_exonic
lalba_exonic_haplo = res_lalba_exonic[[1]]
colnames(lalba_exonic_haplo) <- gsub("group=", "", colnames(lalba_exonic_haplo))
lalba_exonic_haplo[,"Total"] =  round(lalba_exonic_haplo[,"Total"], 2)
lalba_exonic_haplo[,"Desert"] =  round(lalba_exonic_haplo[,"Desert"], 2)
lalba_exonic_haplo[,"Taggar"] =  round(lalba_exonic_haplo[,"Taggar"], 2)
lalba_exonic_haplo[,"Nubian"] =  round(lalba_exonic_haplo[,"Nubian"], 2)
lalba_exonic_haplo[,"Nilotic"] =  round(lalba_exonic_haplo[,"Nilotic"], 2)

notPresent = which(apply(lalba_exonic_haplo[,c("Total", as.character(unique(group)))],1,sum,na.rm=T) == 0)
lalba_exonic_haplo = lalba_exonic_haplo[-notPresent,]
lalba_exonic_haplo = lalba_exonic_haplo[order(lalba_exonic_haplo[, "Total"], decreasing = TRUE),]
rownames(lalba_exonic_haplo) <- paste("Haplotype", LETTERS[1:nrow(lalba_exonic_haplo)])

write.table(lalba_exonic_haplo, "lalba_haplotypes_exonic.txt", sep = "\t", quote=FALSE, na="")

#lgb_upstream
lgb_upstream_haplo = res_lgb_upstream[[1]]
colnames(lgb_upstream_haplo) <- gsub("group=", "", colnames(lgb_upstream_haplo))
lgb_upstream_haplo[,"Total"] =  round(lgb_upstream_haplo[,"Total"], 3)
lgb_upstream_haplo[,"Desert"] =  round(lgb_upstream_haplo[,"Desert"], 3)
lgb_upstream_haplo[,"Taggar"] =  round(lgb_upstream_haplo[,"Taggar"], 3)
lgb_upstream_haplo[,"Nubian"] =  round(lgb_upstream_haplo[,"Nubian"], 3)
lgb_upstream_haplo[,"Nilotic"] =  round(lgb_upstream_haplo[,"Nilotic"], 3)

notPresent = which(apply(lgb_upstream_haplo[,c("Total", as.character(unique(group)))],1,sum,na.rm=T) == 0)
lgb_upstream_haplo = lgb_upstream_haplo[-notPresent,]
lgb_upstream_haplo = lgb_upstream_haplo[order(lgb_upstream_haplo[, "Total"], decreasing = TRUE),]
rownames(lgb_upstream_haplo) <- paste("Haplotype", LETTERS[1:nrow(lgb_upstream_haplo)])

write.table(lgb_upstream_haplo, "lgb_haplotypes_upstream.txt", sep = "\t", quote=FALSE, na="")

#lgb_exonic
lgb_exonic_haplo = res_lgb_exonic[[1]]
colnames(lgb_exonic_haplo) <- gsub("group=", "", colnames(lgb_exonic_haplo))
lgb_exonic_haplo[,"Total"] =  round(lgb_exonic_haplo[,"Total"], 2)
lgb_exonic_haplo[,"Desert"] =  round(lgb_exonic_haplo[,"Desert"], 2)
lgb_exonic_haplo[,"Taggar"] =  round(lgb_exonic_haplo[,"Taggar"], 2)
lgb_exonic_haplo[,"Nubian"] =  round(lgb_exonic_haplo[,"Nubian"], 2)
lgb_exonic_haplo[,"Nilotic"] =  round(lgb_exonic_haplo[,"Nilotic"], 2)

notPresent = which(apply(lgb_exonic_haplo[,c("Total", as.character(unique(group)))],1,sum,na.rm=T) == 0)
lgb_exonic_haplo = lgb_exonic_haplo[-notPresent,]
lgb_exonic_haplo = lgb_exonic_haplo[order(lgb_exonic_haplo[, "Total"], decreasing = TRUE),]
rownames(lgb_exonic_haplo) <- paste("Haplotype", LETTERS[1:nrow(lgb_exonic_haplo)])

write.table(lgb_exonic_haplo, "lgb_haplotypes_exonic.txt", sep = "\t", quote=FALSE, na="")


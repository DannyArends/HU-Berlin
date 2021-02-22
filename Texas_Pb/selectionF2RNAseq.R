
setwd("D:/Edrive/Mouse/Texas_Pb")

top1 <- "gJAX00272462"
top18 <- "gUNC29787730"
topX <- "SXX202459414"


map <- read.table("reblasted_map.txt",sep="\t", header=TRUE, row.names=1)
gts <- read.table("genotypes_F2_filtered_ordered.txt",sep="\t", header=TRUE, row.names=1)
phe <- read.table("F2_phenotypes_cleaned_matched.txt",sep="\t", header=TRUE, row.names=1, na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
rownames(phe) <- gsub("-", ".",rownames(phe))
phe <- phe[colnames(gts),]
map <- map[rownames(gts),]

op <- par(mar = c(10,10,4,2))
heatmap(cor(apply(phe[,c(2,8:22, 25:ncol(phe))],2,as.numeric), use="pair"), margins = c(8,8), scale = "none")

phe <- cbind(BLLadj = NA, phe)

BLLadj <- residuals(lm(BLL ~ sex + urineVolume + waterConsumed, data = phe)) + mean(phe[,"BLL"], na.rm=TRUE)
phe[names(BLLadj),"BLLadj"] <- BLLadj

G1 <- phe[colnames(gts)[which(gts[top1,] == "AA" & gts[top18,] == "AA" & gts[topX,] == "AT")],]
G2 <- phe[colnames(gts)[which(gts[top1,] == "AA" & gts[top18,] == "AG" & gts[topX,] == "AT")],]
G3 <- phe[colnames(gts)[which(gts[top1,] == "AA" & gts[top18,] == "GG" & gts[topX,] == "AT")],]

mean(G1[which(G1[, "sex"] == "F"),1:5][,"BLLadj"],na.rm=T)
mean(G2[which(G2[, "sex"] == "F"),1:5][,"BLLadj"],na.rm=T)
mean(G3[which(G3[, "sex"] == "F"),1:5][,"BLLadj"],na.rm=T)
G1[which(G1[, "sex"] == "F"),1:5]
G2[which(G2[, "sex"] == "F"),1:5]
G3[which(G3[, "sex"] == "F"),1:5]

#Selected F2 animals for gene expression
Low(AA,AA,AT):	
F2.1711_460	25	23.0
F2.1711_636	27	22.0
F2.1711_370	46	34.6

Med(AA,AG,AT):	
F2.1711_292	59	47.9
F2.1711_387	53	41.4
F2.1711_235	54	50.0

High(AA,GG,AT):	
F2.1711_316	100	93.8
F2.1711_664	98	91.0
F2.1711_433	84	78.2


setwd("D:/Edrive/Horse/DNA/")
markerinfo <- read.csv("combined/input/map.txt", sep="\t",colClasses="character")
markerinfo[markerinfo[,"Chr"] == "X","Chr"] <- 32

markerinfo <- cbind(markerinfo, ChrNum = as.numeric(markerinfo[,"Chr"]))
markerinfo <- markerinfo[with(markerinfo, order(ChrNum, MapInfo)), ]


map <- cbind(as.character(markerinfo[, "ChrNum"]), rownames(markerinfo), markerinfo[, "MapInfo"])
write.table(map, "horses.map",sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)

phenotypes <- read.csv("combined/input/phenotypes.txt", sep="\t")
ind <- rownames(phenotypes)[which(!is.na(phenotypes[,"WH"]))]

mm <- which(apply(apply(phenotypes[ind,],1,is.na),1,sum) == 48)
phenotypes <- phenotypes[ind, -mm]

sex     <- as.factor(unlist(phenotypes[,"Sex"]))
strain  <- as.factor(unlist(phenotypes[,"Strain"]))
ageAtM  <- as.numeric(unlist(phenotypes[,"D.Measure"])) - as.numeric(unlist(phenotypes[,"D.Birth"]))
ageAtR  <- as.numeric(unlist(phenotypes[,"D.Racing"])) - as.numeric(unlist(phenotypes[,"D.Birth"]))

phenonames <- c("Age", "WH", "CW", "CH", "NG", "TG", "ChG", "ChD", "ChW", "BLL", "BL", "FCL", "HCL")
racepheno <- c("Distance..km.", "Speed.km.hr.")

pheno <- cbind(FID = rep("Fam0", 48), IID = ind, phenotypes[,phenonames], phenotypes[,racepheno])
write.table(pheno, "horses.pheno", sep=" ", quote=FALSE, row.names=FALSE)

# Calculate the phenotypes after adjusting the fixed effects
phenoC <- matrix(NA, length(c(phenonames, racepheno)), nrow(phenotypes[,]), dimnames = list(c(phenonames, racepheno), rownames(phenotypes)))
for(ph in phenonames) {
  model <- lm(as.numeric(phenotypes[,ph]) ~ sex + ageAtM + strain)
  pCorrected <- model$residuals + model$coefficients["(Intercept)"]
  phenoC[ph, as.numeric(names(pCorrected))] <- round(pCorrected, d=1)
}
for(ph in racepheno) {
  model <- lm(as.numeric(phenotypes[,ph]) ~ sex + ageAtR + strain)
  pCorrected <- model$residuals + model$coefficients["(Intercept)"]
  phenoC[ph, as.numeric(names(pCorrected))] <- round(pCorrected, d=1)
}
phenoCor <- cbind(FID = rep("Fam0", 48), IID = ind, t(phenoC)[,phenonames], t(phenoC)[,racepheno])
write.table(phenoCor, "horses.adjusted.pheno", sep=" ", quote=FALSE, row.names=FALSE)

genotypes_snp <- read.csv("combined/input/genotypes_snp.txt", sep="\t", colClasses="character")
genotypes_snp <- genotypes_snp[rownames(markerinfo), ind]

ped <- cbind(rep("Fam0", 48), ind, rep("0", 48), rep("0", 48), phenotypes[ind,"Sex"])

for(x in 1:nrow(genotypes_snp)){
  a <- unlist(lapply(strsplit(as.character(genotypes_snp[x,]), ""),"[", 1))
  b <- unlist(lapply(strsplit(as.character(genotypes_snp[x,]), ""),"[", 2))
  ped <- cbind(ped, a, b)
}
write.table(ped, "horses.ped",sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE, na="0")

"plink --ped horses.ped --no-pheno --map horses.map --map3 --pheno horses.pheno --mpheno 15"
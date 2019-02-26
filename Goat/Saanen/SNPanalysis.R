### Analysis of A single SNP association in Saanen goats

setwd("D:/Edrive/Goat/DNA/Saanen Association")
genotype <- read.table("Genotypes.txt",sep='\t')
genotype <- genotype[, -c(1,4,5,7,8,9,10,11)]
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]

leistung <- read.table("Leistung.txt", sep = "\t",header=TRUE)
rownames(leistung) <- leistung[,3]

zuchtwerte <- read.table("Zuchtwerte.txt", sep = "\t",header=TRUE)
rownames(zuchtwerte) <- zuchtwerte[,3]

## Check if the owners is associated with the phenotypes
for(phe in zuchtwerte[rownames(genotype),c("ZW.Milch", "ZW.Fett", "ZW.Eiweiss", "B.")]){
  cat(anova(lm(phe ~ zuchtwerte[rownames(genotype),"Besitzer"]))[[5]][1], "\n")
}

## Check if the owners is associated with the phenotypes
for(phe in leistung[rownames(genotype),c("Milch.kg", "Fett..", "Eiw..", "Lakt..", "Milch.kg.1", "Fett...1", "Eiw...1", "Lakt...1")]){
  cat(anova(lm(phe ~ leistung[rownames(genotype),"Besitzer"]))[[5]][1], "\n")
}

owner <- zuchtwerte[rownames(genotype), "Besitzer"]
gt <- factor(genotype[,2], levels=c("AA", "AG", "GG"))
res <- apply(zuchtwerte[rownames(genotype),c("ZW.Milch", "ZW.Fett", "ZW.Eiweiss", "B.")], 2, function(pheno){
  model <- lm(pheno ~ owner + gt)
  return(c(model, anova(model)))
})
lapply(res, "[", 1)

owner <- leistung[rownames(genotype), "Besitzer"]
gt <- factor(genotype[,2], levels=c("AA", "AG", "GG"))
res <- apply(leistung[rownames(genotype),c("Milch.kg", "Fett..", "Eiw..", "Lakt..", "Milch.kg.1", "Fett...1", "Eiw...1", "Lakt...1")], 2, function(pheno){
  model <- lm(pheno ~ owner + gt)
  return(c(model, anova(model)))
})
lapply(res, "[", 1)
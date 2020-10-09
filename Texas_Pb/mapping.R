setwd("D:/Edrive/Mouse/Texas_Pb")

map <- read.table("map_ordered.txt",sep="\t", header=TRUE, row.names=1)
gts <- read.table("genotypes_F2_filtered_ordered.txt",sep="\t", header=TRUE, row.names=1)
phe <- read.table("input/F2_phenotypes.txt",sep="\t", header=TRUE, row.names=1, na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
rownames(phe) <- gsub("-", ".",rownames(phe))
phe <- phe[colnames(gts),]

# Some info

dim(map)
dim(gts)
dim(phe)

BLL <- phe[, "BLL_.ug.mL."]
sex <- phe[, "Sex"]
urineVolume <- phe[, "Urine_Volume_uL"]
ageAtDosing <- phe[, "Age_At_Dosing"]
waterConsumed <- phe[, "Water_Consumed_in_UC_mL"]
weightBeforeDiet <- phe[, "Weight_1_Before.Diet"]

anova(lm(BLL ~ sex)) # highly significant include
anova(lm(BLL ~ sex + urineVolume)) # highly significant include
anova(lm(BLL ~ sex + ageAtDosing)) # significant (might be included)
anova(lm(BLL ~ sex + waterConsumed)) # highly significant include
anova(lm(BLL ~ sex + weightBeforeDiet)) # significant (might be included)

anova(lm(BLL ~ sex + urineVolume + waterConsumed)) # Minimal model for mapping

models <- apply(gts,1,function(gt){
  return(lm(BLL ~ sex + urineVolume + waterConsumed + gt))
})

anovas <- lapply(models, anova)
pvals <- unlist(lapply(anovas, function(x){ return(x["gt", "Pr(>F)"]) }))
lods <- -log10(pvals)
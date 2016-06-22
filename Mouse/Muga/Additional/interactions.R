#
# New QTL mapping of BFMI, using a slightly different population structure correction
#

setwd("E:/Mouse/DNA/MegaMuga")

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genoF2 <- read.table("extra/genotypesF2.txt", sep="\t", check.names=FALSE, colClasses="character")                  # Normal A, H, B genotypes
genoF1 <- read.table("extra/genotypesF1.txt", sep="\t", check.names=FALSE, colClasses="character")                  # Normal A, H, B genotypes
phenotypes <- read.table("extra/phenotypes.txt", sep="\t", check.names=FALSE, colClasses="character")               # Phenotypes

F2  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                       # The F2 individuals
F1  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]                                                      # The F1 individuals

d63 <- as.numeric(phenotypes[F2,"d63"])

res <- matrix(NA, 21260, 21260)

m1i <- 1
apply(genoF2[1:21259, F2],1,function(m1){
  m2i <- m1i + 1
  for(m2i in m1i:21260){
    model <- lm(d63 ~ m1 * genoF2[m2i, F2])
    aovm <- anova(model)[[5]]
    if(length(aovm) == 4) res[m1i,m2i] <- aovm[3]
  })
  m1i <<- m1i + 1
})

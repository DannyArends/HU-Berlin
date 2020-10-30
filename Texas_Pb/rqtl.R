
sex <- phe[,"sex"]
sex[sex == "M"] <- 1
sex[sex == "F"] <- 0

map2 <- map
map2[,2] <- (map2[,2] / 1000000)

# Create a cross object, and see if we can use R/qtl to fix the marker ordering
fake.phe <- round(runif(ncol(gts.num)),1)
cat(paste0(c("fake", "", "", paste0(fake.phe,collapse="\t")), collapse="\t"), "\n", file = "my.cross.csv")
cat(paste0(c("sex", "", "", paste0(sex,collapse="\t")), collapse="\t"), "\n", file = "my.cross.csv",append=TRUE)
cat(paste0(c("pgm", "", "", paste0(rep(1, ncol(gts.num)),collapse="\t")), collapse="\t"), "\n", file = "my.cross.csv",append=TRUE)
write.table(cbind(map2,gts.num), "my.cross.csv",sep="\t", append=TRUE, col.names=FALSE, quote=FALSE)
library(qtl)
mcross <- read.cross(file="my.cross.csv", "csvr", sep = "\t", genotypes=c(-1,0,1))



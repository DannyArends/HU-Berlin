MDCfamily <- "7984"

setwd("D:/Edrive/Mouse/MDC/Genotypes May20")
geno <- read.csv(paste0(MDCfamily, ".txt"), sep = "\t", na.strings = c("", "NA", "-", "X"))

setwd("D:/Edrive/Mouse/MDC/Phenotypes May20")
pheno <- read.csv(paste0(MDCfamily, "_mri.txt"), sep = "\t", na.strings = c("", "NA", "-", "X"))

colnames(geno) <- c("MouseID", "GT", "Uncertain")

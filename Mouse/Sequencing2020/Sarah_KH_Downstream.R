
(128250000 - 127500000) / 242
setwd("D:/Edrive/Mouse/Kurzhaar")

snps <- read.table("KH_Sarah.snps-filtered.vcf")
colnames(snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "KH", "B6-3", "B6-4", "B6-5")
vep <- read.table("KH_Sarah_VEP.vcf")
colnames(vep) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "KH", "B6-3", "B6-4", "B6-5")

info <- unlist(lapply(strsplit(as.character(vep[grep("missense_variant", vep[, "INFO"]),"INFO"]), ";"), "[", 15))
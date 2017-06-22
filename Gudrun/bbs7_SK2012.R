
setwd("D:/")
vcf <- read.table("bbs7.vcf.vep.txt")

bbs7only <- vcf[which(vcf[,6] == "Bbs7"),]


noninteresting <- c("intron_variant", "upstream_gene_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "non_coding_transcript_exon_variant,non_coding_transcript_variant", "intron_variant,non_coding_transcript_variant")

bbs7only[!(bbs7only[,4] %in% noninteresting),-c(6,7,8,10,11,12,13,14,15,16,17,21,22,24,25,27,29:35)]


chr <- "chr3"
start <- 36573142
stop <- 36613477

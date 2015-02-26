# Analysis of Kabadiner horse SNP and performance data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", "ecaballus_gene_ensembl")                                                                    # Biomart for Equus caballus genes

setwd("E:/Horse/DNA/Kabadiner/")
chromosomes <- as.character(c(1:31, "X", "Y", "MT"))
chrinfo     <- read.table("info/chrinfo.txt", sep="\t", header=TRUE)
genotypes   <- read.table("input/cleaned_genotypes.txt", sep = "\t")
map         <- read.table("input/cleaned_map.txt", sep = "\t")

he <- apply(genotypes, 1, function(x){ sum(x == "AG" | x == "AC", na.rm=TRUE) })                                                    # Calculate the number of heterozygous alleles
ho <- apply(genotypes, 1, function(x){ sum(!(x == "AG" | x == "AC"), na.rm=TRUE) })                                                 # Calculate the number of homozygous alleles

window.size   <- 750000
step.size     <- 750000/2
nperms        <- 1000

bins <- NULL
for(chr in chromosomes){
  chr.length <- chrinfo[chrinfo[,"Chr"] == chr, "Length"]                                                                           # Length of the chromosome
  chr.bins <- cbind(Start = seq(1, chr.length, step.size), Stop = seq(1, chr.length, step.size) + window.size)                      # Create our bins
  smap <- map[map[,"Chr"] == chr, ]                                                                                                 # Create a copy of the map (only this chromosome)
  
  stats <- t(apply(chr.bins, 1, function(x){
    chr.markers <- rownames(smap)[as.numeric(smap[,"Position"]) > x["Start"] & as.numeric(smap[,"Position"]) < x["Stop"]]           # Which markers are in my bin
    return(c(length(chr.markers), sum(he[chr.markers]) / sum(ho[chr.markers])))                                                     # Calculate HE / HO ratio
  }))
  chr.bins <- cbind(chr.bins, nSNPs = stats[, 1], Score = stats[, 2])                                                               # Remember the scores

  for(perm in 1:nperms){                                                                                                            # permute this chromosome
    permstats <- t(apply(chr.bins, 1, function(x){
      if(as.numeric(x["nSNPs"]) != 0){                                                                                              # If there are SNPs in the bin
        chr.markers <- rownames(smap)[sample(length(rownames(smap)), as.numeric(x["nSNPs"]))]                                       # Draw N SNPs at random from this chromosome
        return(c(length(chr.markers), sum(he[chr.markers]) / sum(ho[chr.markers])))                                                 # Calculate HE / HO ratio
      }else{
        return(c(0, NaN))                                                                                                           # No SNPs, just return no score
      }
    }))
    chr.bins <- cbind(chr.bins, permstats[, 2])                                                                                     # Add the scores observed during permutation
  }
  cat(paste0("Done with chromosome ",chr,"/",length(chromosomes),"\n"))
  chr.bins <- cbind(Chr = chr, chr.bins)
  bins <- rbind(bins, chr.bins)
}
colnames(bins)[which(colnames(bins) == "")] <- paste0("P" , 1:length(which(colnames(bins) == "")))
write.table(bins,file="analysis/MissingHetro.txt", sep = "\t", row.names=FALSE)

permutationmatrix <- apply(bins[,grep("^P", colnames(bins))],2,as.numeric)

bin.thresholds <- t(apply(permutationmatrix,1,function(x){
  return(c(mean(x) + 5 * sd(x), mean(x) - 5 * sd(x)))
}))

ratio.above <- which(bins[,"Score"] > bin.thresholds[,1])
ratio.below <- which(bins[,"Score"] < bin.thresholds[,2])

regions.above <- apply(bins[ratio.above,1:3], 1, paste0, collapse=":")

regions.below <- apply(bins[ratio.below,1:3], 1, paste0, collapse=":")

bm.above <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene", "external_gene_name", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region", "biotype"), values = list(regions.above,"protein_coding"), mart = bio.mart)
bm.below <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene", "external_gene_name", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region", "biotype"), values = list(regions.below,"protein_coding"), mart = bio.mart)

cat("Found",length(unique(bm.above[,"ensembl_gene_id"])), "genes in hyper-heterozygous regions\n")
cat("Found",length(unique(bm.below[,"ensembl_gene_id"])), "genes in hypo-heterozygous regions\n")

### InnateDB does not have horse GO or Pathway information, we use the human homologs

cat(unique(bm.above[,"hsapiens_homolog_ensembl_gene"]),sep="\n")
# Glycerophospholipid metabolism	3%	0.00408	0.03057
# Tryptophan metabolism	5%	0.00868	0.01953
# Vibrio cholerae infection	4%	0.0149	0.02853
# Glycerolipid metabolism	4%	0.01599	0.02936
# Phosphatidylinositol signaling system	2%	0.0337	0.05616
# Progesterone-mediated oocyte maturation	2%	0.04075	0.06548

cat(unique(bm.below[,"hsapiens_homolog_ensembl_gene"]),sep="\n")




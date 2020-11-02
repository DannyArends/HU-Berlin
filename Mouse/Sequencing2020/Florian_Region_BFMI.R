# Call SNPs for Kourosh using all mouse strains

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}

callSNPs <- function(bamfiles, chr = 1, startpos = 1, endpos = 2, outname = "mySNPs") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome
  region = paste0(chr, ":", format(startpos, scientific = FALSE), "-", format(endpos, scientific = FALSE)) # Region requested

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ", region, " -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && DP>10' - -o ", outname, ".snps-filtered.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

regions <- rbind(c(3,36599339,36599423),
      c(3,36602829,36602959),
      c(5,142905564,142905692),
      c(5,142906652,142906754),
      c(6,125162025,125165291),
      c(6,125165271,125165311),
      c(6,90557202,90557351),
      c(6,90559242,90559476),
      c(8,123420607,123422010),
      c(18,82572825,82572926),
      c(18,82575172,82575207),
      c(18,82578936,82579058))

for(r in 1:nrow(regions)){
callSNPs(c("/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/10/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/11/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/12/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/13/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/14/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam"
           ), regions[r,1], regions[r,2], regions[r,3], paste0("bfmi:", regions[r,1], ":",regions[r,2],"-", regions[r,3]))
}

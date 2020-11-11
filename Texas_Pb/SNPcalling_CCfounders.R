# Call SNPs for Texas_Pb using all CC founder strains

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

founders <- c("/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/A_J.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129S1_SvImJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NOD_ShiLtJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZO_HlLtJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/CAST_EiJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/PWK_PhJ.bam",
              "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/WSB_EiJ.bam")

callSNPs(founders, 1, 124548738, 184926264, "CCfounders_QTL1")
callSNPs(founders, 1, 135464050, 157672910, "CCfounders_QTL2")
callSNPs(founders, 3, 85146443, 147514132, "CCfounders_QTL3")
callSNPs(founders, 4, 58831312, 88222686, "CCfounders_QTL4")
callSNPs(founders, 6, 76940963, 114257130, "CCfounders_QTL5")
callSNPs(founders, 7, 62305639, 81923457, "CCfounders_QTL6")
callSNPs(founders, 7, 62305639, 97260868, "CCfounders_QTL7")
callSNPs(founders, 7, 19482853, 119720011, "CCfounders_QTL8")
callSNPs(founders, 8, 31799796, 104944836, "CCfounders_QTL9")

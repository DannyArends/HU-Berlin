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

samples <- c("/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/1/SJLP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/2/NZOP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/3/NZOP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/10/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/11/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/12/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/13/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/14/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/15/B6-3P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/16/B6-4P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/17/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/18/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/19/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/20/AKRP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/21/AKRP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/22/AKRP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/AKR_J.bam",
           "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
           "/home/danny/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZO_HlLtJ.bam" )

callSNPs(samples, 1,  74938816, 74940881, "mmu-miR-375-5p")
callSNPs(samples, 15, 73125663,73127732, "mmu-miR-151-3p")
callSNPs(samples, 12, 36865203, 36867278, "mmu-miR-5099")
callSNPs(samples, 4,  36667509, 36669587, "mmu-miR-873a-3p")


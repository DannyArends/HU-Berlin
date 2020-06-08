#KH mouse chr 7: 127 500 000 to 128 250 000

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
  bcftools = "/home/arends/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/arends/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome
  region = paste0(chr, ":", format(startpos, scientific = FALSE), "-", format(endpos, scientific = FALSE)) # Region requested

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ", region, " -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -v snps -i '%QUAL>=30 && DP>10' - -o ", outname, ".snps-filtered.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

callSNPs(c("/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/4/KHP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/5/KHP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/15/B6-3P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/16/B6-4P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/17/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/18/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/19/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam"
           ), 7, 127500000, 128250000, "KH_Sarah")


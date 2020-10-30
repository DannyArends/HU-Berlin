# Call SNPs for Florian using BFMI

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
  cmd3 <- paste0(bcftools, " view -i '%QUAL>=30 && DP>10' - -o ", outname, ".snps-filtered.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

bams <- c("/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/10/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
          "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/11/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
          "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/12/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
          "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/13/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
          "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/14/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam")

callSNPs(bams, 6, 90559981, 90559991, "Aldh1l1_1")
callSNPs(bams, 6, 90562619, 90562627, "Aldh1l1_2")
callSNPs(bams, 6, 90563396, 90563415, "Aldh1l1_3")

callSNPs(bams, 11, 102896857, 102896886, "Gfap_1")
callSNPs(bams, 11, 102896666, 102896670, "Gfap_2")
callSNPs(bams, 11, 102895794, 102895810, "Gfap_3")

callSNPs(bams, 16, 91226440, 91226459, "Olig2_1")
callSNPs(bams, 16, 91226626, 91226645, "Olig2_2")

callSNPs(bams, 11, 118909304, 118909312, "Ribfox3_1")
callSNPs(bams, 11, 118811406, 118811416, "Ribfox3_2")
callSNPs(bams, 11, 118811278, 118811287, "Ribfox3_3")
callSNPs(bams, 11, 118670744, 118670753, "Ribfox3_4")

callSNPs(bams, 8, 89044030, 89044037, "Sall1_1")
callSNPs(bams, 8, 89042433, 89042444, "Sall1_2")
callSNPs(bams, 8, 89033379, 89033398, "Sall1_3")

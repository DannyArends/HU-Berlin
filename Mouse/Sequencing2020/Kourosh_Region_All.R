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
  bcftools = "/home/arends/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/arends/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa" #Reference genome
  region = paste0(chr, ":", format(startpos, scientific = FALSE), "-", format(endpos, scientific = FALSE)) # Region requested

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ", region, " -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -v snps -i '%QUAL>=30 && DP>10' - -o ", outname, ".snps-filtered.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

callSNPs(c("/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/1/SJLP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/2/NZOP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/3/NZOP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/4/KHP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/5/KHP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/6/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/7/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/8/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/9/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/10/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/11/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/12/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/13/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/14/BFMI860-12P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/15/B6-3P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/16/B6-4P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/17/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/18/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/19/B6-5P_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/20/AKRP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/21/AKRP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/22/AKRP_trimmed.aligned.sorted.dedup.recalibrated.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129P2_OlaHsd.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129S1_SvImJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129S5SvEvBrd.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/A_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/AKR_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/BALB_cJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/BTBR_T__Itpr3tf_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/BUB_BnJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeH.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_10J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BR_cdJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57L_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C58_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/CAST_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/CBA_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_1J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/FVB_NJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/I_LnJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/JF1_MsJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/KK_HiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/LEWES_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/LG_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/LP_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/MOLF_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NOD_ShiLtJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZB_B1NJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZO_HlLtJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZW_LacJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/PWK_PhJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/RF_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SEA_GnJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SJL_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SM_J.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SPRET_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/ST_bJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/WSB_EiJ.bam",
           "/home/arends/NAS/old_naszb/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/ZALENDE_EiJ.bam"
           ), 3, 36000000, 37000000, "MGP_Chr3:36000000-37000000")


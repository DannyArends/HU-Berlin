setwd("/home/arends/MouseAlignment2020/")

allfiles <- rbind(
# SJL
 c("SJL", "/home/arends/NAS/Mouse/DNA/Sequencing/SJL_Dife_2016/ext_L7255-3_SJL_S16_R1_001.fastq.gz", "/home/arends/NAS/Mouse/DNA/Sequencing/SJL_Dife_2016/ext_L7255-3_SJL_S16_R2_001.fastq.gz"),
# NZO
 c("NZO", "/home/arends/NAS/Mouse/DNA/Sequencing/NZO_Dife_2016/ext_L7254-1_NZO_S1_R1_001.fastq.gz", "/home/arends/NAS/Mouse/DNA/Sequencing/NZO_Dife_2016/ext_L7254-1_NZO_S1_R2_001.fastq.gz"),
 c("NZO", "/home/arends/NAS/Mouse/DNA/Sequencing/NZO_Dife_2016/ext_L7254-2_NZO_S3_R1_001.fastq.gz", "/home/arends/NAS/Mouse/DNA/Sequencing/NZO_Dife_2016/ext_L7254-2_NZO_S3_R2_001.fastq.gz"),
# KH Mouse
 c("KH", "/home/arends/NAS/Mouse/DNA/Sequencing/KH_2020/fastq/M_KH_2_FDSW202341665-1r_H35JYDSXY_L1_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/KH_2020/fastq/M_KH_2_FDSW202341665-1r_H35JYDSXY_L1_2.fq.gz"),
 c("KH", "/home/arends/NAS/Mouse/DNA/Sequencing/KH_2020/fastq/M_KH_2_FDSW202341665-1r_H35KFDSXY_L3_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/KH_2020/fastq/M_KH_2_FDSW202341665-1r_H35KFDSXY_L3_2.fq.gz"),
# BFMI861-S1
 c("BFMI861-S1", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S1_2016/ext_L7256-1_BFMI861-S1_S2_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S1_2016/ext_L7256-1_BFMI861-S1_S2_R2_001.fastq.gz"),
 c("BFMI861-S1", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S1_2016/ext_L7256-2_BFMI861-S1_S4_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S1_2016/ext_L7256-2_BFMI861-S1_S4_R2_001.fastq.gz"),
# BFMI-861-S2
 c("BFMI861-S2", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S2_2016/ext_L7257-1_BFMI861-S2_S5_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S2_2016/ext_L7257-1_BFMI861-S2_S5_R2_001.fastq.gz"),
 c("BFMI861-S2", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S2_2016/ext_L7257-2_BFMI861-S2_S6_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI861-S2_2016/ext_L7257-2_BFMI861-S2_S6_R2_001.fastq.gz"),
# BFMI860-S12
 c("BFMI860-12-2014", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860_MDC_2014/BFMI860_L008_R1.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860_MDC_2014/BFMI860_L008_R2.fastq.gz"),
 c("BFMI860-12-2016", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2016/ext_L7258-2_BFMI860-S12_S2_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2016/ext_L7258-2_BFMI860-S12_S2_R2_001.fastq.gz"),
 c("BFMI860-12-2016", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2016/ext_L7258-1_BFMI860-S12_S1_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2016/ext_L7258-1_BFMI860-S12_S1_R2_001.fastq.gz"),
 c("BFMI860-12-2016", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2019/ext_L7258-1_BFMI860-S12_S9_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2019/ext_L7258-1_BFMI860-S12_S9_R2_001.fastq.gz"),
 c("BFMI860-12-2016", "/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2019/ext_L7258-1_BFMI860-S12_S23_R1_001.fastq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/BFMI860-S12_2019/ext_L7258-1_BFMI860-S12_S23_R2_001.fastq.gz"),
# B6
 c("B6-3", "/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_3_FDSW202341661-1r_H35G2DSXY_L2_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_3_FDSW202341661-1r_H35G2DSXY_L2_2.fq.gz"),
 c("B6-4", "/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_4_FDSW202341662-1r_H35G2DSXY_L2_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_4_FDSW202341662-1r_H35G2DSXY_L2_2.fq.gz"),
 c("B6-5", "/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35JYDSXY_L1_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35JYDSXY_L1_2.fq.gz"),
 c("B6-5", "/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35KFDSXY_L3_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35KFDSXY_L3_2.fq.gz"),
 c("B6-5", "/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35MNDSXY_L1_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35MNDSXY_L1_2.fq.gz"),
# AKR
 c("AKR", "/home/arends/NAS/Mouse/DNA/Sequencing/AKR_2020/fastq/M_AKR_1_FDSW202341666-1r_H35JYDSXY_L1_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/AKR_2020/fastq/M_AKR_1_FDSW202341666-1r_H35JYDSXY_L1_2.fq.gz"),
 c("AKR", "/home/arends/NAS/Mouse/DNA/Sequencing/AKR_2020/fastq/M_AKR_1_FDSW202341666-1r_H35KFDSXY_L3_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/AKR_2020/fastq/M_AKR_1_FDSW202341666-1r_H35KFDSXY_L3_2.fq.gz"),
 c("AKR", "/home/arends/NAS/Mouse/DNA/Sequencing/AKR_2020/fastq/M_AKR_1_FDSW202341666-1r_H35MNDSXY_L1_1.fq.gz","/home/arends/NAS/Mouse/DNA/Sequencing/AKR_2020/fastq/M_AKR_1_FDSW202341666-1r_H35MNDSXY_L1_2.fq.gz")
)

colnames(allfiles) <- c("Sample", "R1", "R2")
write.table(allfiles, "inputfiles.txt", sep = "\t", quote=FALSE, row.names = FALSE)

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

for(x in 1:nrow(allfiles)){
  command <- paste0("Rscript pipeline.R ", x, " ", allfiles[x,1], " ", allfiles[x,2], " ", allfiles[x,3])
  execute(command, TRUE)
}



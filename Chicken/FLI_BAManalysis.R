# Get the fasta & DBSNP file
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/galGal4/bigZips/galGal4.fa.gz
wget ftp://ftp.ensembl.org/pub/release-88/variation/vcf/gallus_gallus/Gallus_gallus.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-88/variation/vcf/gallus_gallus/Gallus_gallus.vcf.gz.tbi
#unzip the reference
gunzip galGal4.fa.gz

#index the reference
samtools faidx galGal4.fa
java -jar /home/shijie/tools/picard-tools-1.119/CreateSequenceDictionary.jar R=galGal4.fa O=galGal4.dict


floc <- "/home/shijie/Data_FLI/DNA_Seq_BAM/"

setwd(floc)

files <- dir()[grep("bam", dir())]
indexes <- gsub("bam", "bai", files)
files <- paste0(floc, files)

setwd("~/ChickenDNASeq")
# Index the bam files
for(x in 1:length(files)) {
  system(paste0("samtools index ", files[x], " ", indexes[x]), intern=FALSE)
}

cat(paste0("-I ", files, collapse=" \\\n "), file="tmp.txt")  # Create the long list of input files needed for GATK


nohup java -Xmx12g -jar ~/tools/GATK/GenomeAnalysisTK.jar -nct 28 -T HaplotypeCaller \
-R ~/ChickenDNASeq/reference/galGal4.fa \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_AB_0001_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_AR_0002_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_AS_0003_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_BA_0004_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_BH_0006_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_BK_0005_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_CG_0007_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_CS_0008_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_CW_0037_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_DG_0009_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_DL_0038_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_FG_0010_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_HO_0011_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_IT_0015_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_KA_0016_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_KG_0017_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_KS_0039_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_KW_0014_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_LE_0019_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_MA_0020_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_MR_0021_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_NH_0025_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_OF_0040_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_OH_0022_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_OM_0024_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_OR_0023_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_PA_0018_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_SA_0029_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_SB_0026_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_SE_0027_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_SH_0028_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_SN_0030_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_TO_0031_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_WT_0032_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_WY_0033_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_YO_0034_4_68600000_82000000.mdup.sorted.merged.bam \
-I /home/shijie/Data_FLI/DNA_Seq_BAM/pl_ZC_0035_4_68600000_82000000.mdup.sorted.merged.bam \
--dbsnp ~/ChickenDNASeq/reference/Gallus_gallus.vcf.gz \
-L targets.interval_list \
-dontUseSoftClippedBases \
-o SNPsFLIDNASeq.vcf &

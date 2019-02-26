
## Dump the aligned reads of SRR5266611 in a region of the genome (3:36573142-36613477) to a BAM file

getBAM <- function(sraID, out = "~/dump/sam/", region = NULL, header = TRUE, minMapQ = 0, reconstruct = FALSE){
  cmd = "./sam-dump"
  outloc = paste0(out, sraID, ".bam")
  if(!is.null(region)) cmd = paste(cmd, "--aligned-region", region)
  if(!header) cmd = paste(cmd, "--no-header")
  if(reconstruct) cmd = paste(cmd, "--header")
  if(minMapQ > 0) cmd = paste(cmd, "--min-mapq", minMapQ)
  cmd = paste(cmd, sraID, "| samtools view -bS - >", outloc)
  return(cmd)
}

## Retina
./sam-dump --aligned-region 3:36272064-36712399 SRX2173048 | samtools view -bS - > ~/fastq-dump/SRX2173048.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173049 | samtools view -bS - > ~/fastq-dump/SRX2173049.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173050 | samtools view -bS - > ~/fastq-dump/SRX2173050.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173051 | samtools view -bS - > ~/fastq-dump/SRX2173051.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173052 | samtools view -bS - > ~/fastq-dump/SRX2173052.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173053 | samtools view -bS - > ~/fastq-dump/SRX2173053.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173054 | samtools view -bS - > ~/fastq-dump/SRX2173054.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173055 | samtools view -bS - > ~/fastq-dump/SRX2173055.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173056 | samtools view -bS - > ~/fastq-dump/SRX2173056.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173057 | samtools view -bS - > ~/fastq-dump/SRX2173057.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173058 | samtools view -bS - > ~/fastq-dump/SRX2173058.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173059 | samtools view -bS - > ~/fastq-dump/SRX2173059.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173060 | samtools view -bS - > ~/fastq-dump/SRX2173060.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173061 | samtools view -bS - > ~/fastq-dump/SRX2173061.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173062 | samtools view -bS - > ~/fastq-dump/SRX2173062.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173063 | samtools view -bS - > ~/fastq-dump/SRX2173063.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173064 | samtools view -bS - > ~/fastq-dump/SRX2173064.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173065 | samtools view -bS - > ~/fastq-dump/SRX2173065.bam
./sam-dump --aligned-region 3:36272064-36712399 SRX2173066 | samtools view -bS - > ~/fastq-dump/SRX2173066.bam

samtools index SRX2173048.bam
samtools index SRX2173049.bam
samtools index SRX2173050.bam
samtools index SRX2173051.bam
samtools index SRX2173052.bam
samtools index SRX2173053.bam
samtools index SRX2173054.bam
samtools index SRX2173055.bam
samtools index SRX2173056.bam
samtools index SRX2173057.bam
samtools index SRX2173058.bam
samtools index SRX2173059.bam
samtools index SRX2173060.bam
samtools index SRX2173061.bam
samtools index SRX2173062.bam
samtools index SRX2173063.bam
samtools index SRX2173064.bam
samtools index SRX2173065.bam
samtools index SRX2173066.bam

## Testis
./sam-dump --aligned-region chr3:36272064-36712399 SRX4004711 | samtools view -bS - > ~/fastq-dump/SRX4004711.bam
./sam-dump --aligned-region chr3:36272064-36712399 SRX4004710 | samtools view -bS - > ~/fastq-dump/SRX4004710.bam
./sam-dump --aligned-region chr3:36272064-36712399 SRX4004713 | samtools view -bS - > ~/fastq-dump/SRX4004713.bam

samtools index SRX4004711.bam
samtools index SRX4004710.bam
samtools index SRX4004713.bam

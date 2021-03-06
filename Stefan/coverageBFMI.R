
seqtk sample -s100 /home/neubert/BFMImm10/BFMI860_L008_R1.fastq.gz 10000 > BFMI_subset_R1.fq
seqtk sample -s100 /home/neubert/BFMImm10/BFMI860_L008_R2.fastq.gz 10000 > BFMI_subset_R2.fq


## Trimming the reads

java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar PE BFMI_subset_R1.fq BFMI_subset_R2.fq BFMI_subset_R1.P.fastq.gz BFMI_subset_R1.U.fastq.gz BFMI_subset_R2.P.fastq.gz BFMI_subset_R2.U.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar PE /home/neubert/BFMImm10/BFMI860_L008_R1.fastq.gz /home/neubert/BFMImm10/BFMI860_L008_R2.fastq.gz BFMI_R1.P.fastq.gz BFMI_R1.U.fastq.gz BFMI_R2.P.fastq.gz BFMI_R2.U.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# Creating dictionaries for the reference
java -jar /opt/picard-tools-1.99/CreateSequenceDictionary.jar R=/home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa O=/home/share/genomes/mm10/Mus_musculus.GRCm38.dna.dict
bowtie2-build -f /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.bt2idx
samtools faidx /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa

tabix /home/share/genomes/mm10/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz
tabix /home/share/genomes/mm10/mgp.v5.merged.snps_all.dbSNP142.vcf.gz 

# Alignment of the reads to the reference

#bowtie2 -x /home/arends/RNASeq/Reference/Mus_musculus.GRCm38.dna.bt2idx -1 BFMI_R1.P.fastq.gz -2 BFMI_R2.P.fastq.gz -U BFMI_R1.U.fastq.gz,BFMI_R2.U.fastq.gz -X 2000 -I 50 -p 4 | samtools view -bS - > BFMI_aligned.bam
bowtie2 -x /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.bt2idx -1 BFMI_subset_R1.P.fastq.gz -2 BFMI_subset_R2.P.fastq.gz -U BFMI_subset_R1.U.fastq.gz,BFMI_subset_R2.U.fastq.gz -X 2000 -I 50 | samtools view -bS - > BFMI_subset_aligned.bam

# Create some flagstats for the alignment
nohup samtools sort -@ 4 -m 2G BFMI_subset_aligned.bam BFMI_aligned_sorted &
java -jar /opt/picard-tools-1.99/MarkDuplicates.jar REMOVE_DUPLICATES=true INPUT=BFMI_aligned_sorted.bam OUTPUT=BFMI_aligned_sorted_rmdup.bam METRICS_FILE=BFMI_metrics.txt
samtools flagstat BFMI_aligned_sorted_rmdup.bam

#Add read groups
java -Xmx4g -jar /opt/picard-tools-1.99/AddOrReplaceReadGroups.jar INPUT=BFMI_aligned_sorted_rmdup.bam OUTPUT=BFMI_aligned_sorted_rmdup_rg.bam CREATE_INDEX=false RGID=860 RGLB=LIB860 RGPL=Illumina RGPU=X RGSM=860
samtools index BFMI_aligned_sorted_rmdup_rg.bam

### indel realignment 
#Only once
java -Xmx4g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa -known /home/share/genomes/mm10/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz -o /home/share/genomes/mm10/output.indels.intervals -U ALLOW_N_CIGAR_READS

# Realign
java -Xmx4g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa -targetIntervals /home/arends/DNASeq/output.indels.intervals -maxReads 150000 -I BFMI_aligned_sorted_rmdup_rg.bam -o BFMI_aligned_sorted_rmdup_rg_realigned.bam -known /home/share/genomes/mm10/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz -U ALLOW_N_CIGAR_READS

### Base Recalibration
java -Xmx4g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -nct 4 -T BaseRecalibrator -R /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa -knownSites /home/share/genomes/mm10/mgp.v5.merged.snps_all.dbSNP142.vcf.gz -I BFMI_aligned_sorted_rmdup_rg_realigned.bam -o BFMI_1.covariates -U ALLOW_N_CIGAR_READS
java -Xmx4g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -nct 4 -T PrintReads -R /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa -I BFMI_aligned_sorted_rmdup_rg_realigned.bam -BQSR BFMI_1.covariates -U ALLOW_N_CIGAR_READS -o BFMI_aligned_sorted_rmdup_rg_realigned_recalibrated.bam
java -Xmx4g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -nct 4 -T BaseRecalibrator -R /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa -knownSites /home/share/genomes/mm10/mgp.v5.merged.snps_all.dbSNP142.vcf.gz  -I BFMI_aligned_sorted_rmdup_rg_realigned_recalibrated.bam -o BFMI_2.covariates -U ALLOW_N_CIGAR_READS
java -Xmx4g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T AnalyzeCovariates -R /home/share/genomes/mm10/Mus_musculus.GRCm38.dna.fa -before BFMI_1.covariates -after BFMI_2.covariates -U ALLOW_N_CIGAR_READS -plots BFMI_recalibration.pdf


mgp.v5.merged.snps_all.dbSNP142.vcf





chr <- 3
pos <- 14846228

command <- paste0("samtools mpileup -f /home/neubert/ref_genomes/mm10/Mus_musculus.GRCm38.dna.fa -r chr",chr,":",pos,"-",pos," /home/neubert/alt_BFMImm10SK2014/Serverdateienzwischenstop/860v2.pes.s.rmd.rg.realigned.recal.bam")
system(command, intern = TRUE)
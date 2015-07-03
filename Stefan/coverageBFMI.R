

seqtk sample -s100 /home/neubert/BFMImm10/BFMI860_L008_R1.fastq.gz 10000 > BFMI_subset_R1.fq
seqtk sample -s100 /home/neubert/BFMImm10/BFMI860_L008_R2.fastq.gz 10000 > BFMI_subset_R2.fq

# java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar PE BFMI_subset_R1.fq BFMI_subset_R2.fq BFMI_subset_R1.P.fastq.gz BFMI_subset_R1.U.fastq.gz BFMI_subset_R2.P.fastq.gz BFMI_subset_R2.U.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar PE /home/neubert/BFMImm10/BFMI860_L008_R1.fastq.gz /home/neubert/BFMImm10/BFMI860_L008_R2.fastq.gz BFMI_R1.P.fastq.gz BFMI_R1.U.fastq.gz BFMI_R2.P.fastq.gz BFMI_R2.U.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


bowtie2-build -f Mus_musculus.GRCm38.dna.fa Mus_musculus.GRCm38.dna.bt2idx

bowtie2 -x /home/arends/RNASeq/Reference/Mus_musculus.GRCm38.dna.bt2idx -1 BFMI_R1.P.fastq.gz -2 BFMI_R2.P.fastq.gz -S BFMIaligned.sam

chr <- 3
pos <- 14846228

command <- paste0("samtools mpileup -f /home/neubert/ref_genomes/mm10/Mus_musculus.GRCm38.dna.fa -r chr",chr,":",pos,"-",pos," /home/neubert/alt_BFMImm10SK2014/Serverdateienzwischenstop/860v2.pes.s.rmd.rg.realigned.recal.bam")
system(command, intern = TRUE)
>mmu-mir-5099 (chr12:36816205-36816278) ENSMUSG00000092998 
ttagatcgatgtggtgctcc 
>mmu-mir-375-5p (chr1:74900658-74900721) ENSMUSG00000065616 
gcgacgagcccctcgcacaaac
>mmu-mir-151-3p (chr15:73254815-73254882) ENSMUSG00000065612 
ctagactgaggctccttgagg
>mmu-mir-873a-3p (chr4:36668510-36668586) ENSMUSG00000078057 
gagactgacaagttcccggga

# Make fasta file from fastq reads
gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-1_860-29698/HF3W7DRXY/L001/BFMI-1_860-29698_S9_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-1_860-29698/HF3W7DRXY/L001/BFMI-1_860-29698_S9_L001_R1_001.fa
gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-2_860-29699/HF3W7DRXY/L001/BFMI-2_860-29699_S10_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-2_860-29699/HF3W7DRXY/L001/BFMI-2_860-29699_S10_L001_R1_001.fa
gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-3_860-29700/HF3W7DRXY/L001/BFMI-3_860-29700_S11_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-3_860-29700/HF3W7DRXY/L001/BFMI-3_860-29700_S11_L001_R1_001.fa
gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-4_860-29692/HF3W7DRXY/L001/BFMI-4_860-29692_S12_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-4_860-29692/HF3W7DRXY/L001/BFMI-4_860-29692_S12_L001_R1_001.fa

gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/B6-1_1009437/HF3W7DRXY/L001/B6-1_1009437_S5_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/B6-1_1009437/HF3W7DRXY/L001/B6-1_1009437_S5_L001_R1_001.fa
gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/B6-2_1009438/HF3W7DRXY/L001/B6-2_1009438_S6_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/B6-2_1009438/HF3W7DRXY/L001/B6-2_1009438_S6_L001_R1_001.fa
gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/B6-3_1009442/HF3W7DRXY/L001/B6-3_1009442_S7_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/B6-3_1009442/HF3W7DRXY/L001/B6-3_1009442_S7_L001_R1_001.fa
gunzip -c /home/danny/NAS/Mouse/RNAseqKourosh/B6-4_1009443/HF3W7DRXY/L001/B6-4_1009443_S8_L001_R1_001.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/danny/NAS/Mouse/RNAseqKourosh/B6-4_1009443/HF3W7DRXY/L001/B6-4_1009443_S8_L001_R1_001.fa

# Get the number of reads per file (without accounting for duplicates)
zcat /home/danny/NAS/Mouse/RNAseqKourosh/B6-1_1009437/HF3W7DRXY/L001/B6-1_1009437_S5_L001_R1_001.fastq.gz |wc -l
zcat /home/danny/NAS/Mouse/RNAseqKourosh/B6-2_1009438/HF3W7DRXY/L001/B6-2_1009438_S6_L001_R1_001.fastq.gz |wc -l
zcat /home/danny/NAS/Mouse/RNAseqKourosh/B6-3_1009442/HF3W7DRXY/L001/B6-3_1009442_S7_L001_R1_001.fastq.gz |wc -l
zcat /home/danny/NAS/Mouse/RNAseqKourosh/B6-4_1009443/HF3W7DRXY/L001/B6-4_1009443_S8_L001_R1_001.fastq.gz |wc -l

zcat /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-1_860-29698/HF3W7DRXY/L001/BFMI-1_860-29698_S9_L001_R1_001.fastq.gz |wc -l
zcat /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-2_860-29699/HF3W7DRXY/L001/BFMI-2_860-29699_S10_L001_R1_001.fastq.gz |wc -l
zcat /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-3_860-29700/HF3W7DRXY/L001/BFMI-3_860-29700_S11_L001_R1_001.fastq.gz |wc -l
zcat /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-4_860-29692/HF3W7DRXY/L001/BFMI-4_860-29692_S12_L001_R1_001.fastq.gz |wc -l

# Make a blastDB out of the reads so we can search for our miRNAs
makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-1_860-29698/HF3W7DRXY/L001/BFMI-1_860-29698_S9_L001_R1_001.fa -dbtype nucl
makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-2_860-29699/HF3W7DRXY/L001/BFMI-2_860-29699_S10_L001_R1_001.fa -dbtype nucl
makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-3_860-29700/HF3W7DRXY/L001/BFMI-3_860-29700_S11_L001_R1_001.fa -dbtype nucl
makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-4_860-29692/HF3W7DRXY/L001/BFMI-4_860-29692_S12_L001_R1_001.fa -dbtype nucl

makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/B6-1_1009437/HF3W7DRXY/L001/B6-1_1009437_S5_L001_R1_001.fa -dbtype nucl
makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/B6-2_1009438/HF3W7DRXY/L001/B6-2_1009438_S6_L001_R1_001.fa -dbtype nucl
makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/B6-3_1009442/HF3W7DRXY/L001/B6-3_1009442_S7_L001_R1_001.fa -dbtype nucl
makeblastdb -in /home/danny/NAS/Mouse/RNAseqKourosh/B6-4_1009443/HF3W7DRXY/L001/B6-4_1009443_S8_L001_R1_001.fa -dbtype nucl

# Blast against our databases
blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-1_860-29698/HF3W7DRXY/L001/BFMI-1_860-29698_S9_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_BFMI860-1.txt
blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-2_860-29699/HF3W7DRXY/L001/BFMI-2_860-29699_S10_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_BFMI860-2.txt
blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-3_860-29700/HF3W7DRXY/L001/BFMI-3_860-29700_S11_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_BFMI860-3.txt
blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/BFMI-4_860-29692/HF3W7DRXY/L001/BFMI-4_860-29692_S12_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_BFMI860-4.txt

blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/B6-1_1009437/HF3W7DRXY/L001/B6-1_1009437_S5_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_B6-1.txt
blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/B6-2_1009438/HF3W7DRXY/L001/B6-2_1009438_S6_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_B6-2.txt
blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/B6-3_1009442/HF3W7DRXY/L001/B6-3_1009442_S7_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_B6-3.txt
blastn -task blastn -query /home/danny/NAS/Mouse/RNAseqKourosh/miRNA.fasta -db /home/danny/NAS/Mouse/RNAseqKourosh/B6-4_1009443/HF3W7DRXY/L001/B6-4_1009443_S8_L001_R1_001.fa -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 100000 -out /home/danny/NAS/Mouse/RNAseqKourosh/blast_miRNA_B6-4.txt


cat Mus_musculus.GRCm38.81.gtf | grep "miRNA" > allMiRNA.gtf

setwd("D:/Edrive/Mouse/miRNA Kourosh")
total <- read.table("total_reads.txt", sep = "\t", row.names=1)

tt <- rbind("B6-1" = c(1000 * table(read.table("blast_miRNA_B6-1.txt", sep = "\t")[,1]) / total["B6-1",2]),
            "B6-2" = c(1000 * table(read.table("blast_miRNA_B6-2.txt", sep = "\t")[,1]) / total["B6-2",2]),
            "B6-3" = c(1000 * table(read.table("blast_miRNA_B6-3.txt", sep = "\t")[,1]) / total["B6-3",2]),
            "B6-4" = c(1000 * table(read.table("blast_miRNA_B6-4.txt", sep = "\t")[,1]) / total["B6-4",2]),
            "BFMI860-1" = c(1000 * table(read.table("blast_miRNA_BFMI860-1.txt", sep = "\t")[,1]) / total["BFMI860-1",2]),
            "BFMI860-2" = c(1000 * table(read.table("blast_miRNA_BFMI860-2.txt", sep = "\t")[,1]) / total["BFMI860-2",2]),
            "BFMI860-3" = c(1000 * table(read.table("blast_miRNA_BFMI860-3.txt", sep = "\t")[,1]) / total["BFMI860-3",2]),
            "BFMI860-4" = c(1000 * table(read.table("blast_miRNA_BFMI860-4.txt", sep = "\t")[,1]) / total["BFMI860-4",2]))

tt

tt1 <- rbind("B6-1" = c(1000 * table(read.table("blast_miRNA_B6-1.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_B6-1.txt", sep = "\t")[,1]))),
            "B6-2" = c(1000 * table(read.table("blast_miRNA_B6-2.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_B6-2.txt", sep = "\t")[,1]))),
            "B6-3" = c(1000 * table(read.table("blast_miRNA_B6-3.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_B6-3.txt", sep = "\t")[,1]))),
            "B6-4" = c(1000 * table(read.table("blast_miRNA_B6-4.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_B6-4.txt", sep = "\t")[,1]))),
            "BFMI860-1" = c(1000 * table(read.table("blast_miRNA_BFMI860-1.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_BFMI860-1.txt", sep = "\t")[,1]))),
            "BFMI860-2" = c(1000 * table(read.table("blast_miRNA_BFMI860-2.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_BFMI860-2.txt", sep = "\t")[,1]))),
            "BFMI860-3" = c(1000 * table(read.table("blast_miRNA_BFMI860-3.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_BFMI860-3.txt", sep = "\t")[,1]))),
            "BFMI860-4" = c(1000 * table(read.table("blast_miRNA_BFMI860-4.txt", sep = "\t")[,1]) / sum(table(read.table("blast_miRNA_BFMI860-4.txt", sep = "\t")[,1]))))

tt1




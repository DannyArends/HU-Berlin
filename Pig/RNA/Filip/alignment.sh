

# Full Genome
wget ftp://ftp.ensembl.org/pub/release-95/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
# Full Transcriptome
wget ftp://ftp.ensembl.org/pub/release-95/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.95.gtf.gz

~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.gtf.gz


cd ~/genomes/SScrofa/
gunzip --keep Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz


# Create some required indexing files
java -Xmx16g -jar /home/danny/picard-tools-1.119/CreateSequenceDictionary.jar R=~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa O=~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.dict
/home/danny/bowtie2-2.3.2/bowtie2-build -f ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx
/home/danny/Github/samtools/samtools faidx ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

# Link the fasta for bowtie2
ln -s Sus_scrofa.Sscrofa11.1.dna.toplevel.fa Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx.fa

# Create the transcriptome index file
/home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -G Sus_scrofa.Sscrofa11.1.95.gtf --transcriptome-index=Sus_scrofa.Sscrofa11.1.95.genes ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx


cd ~/NAS/Pig/S906\ wholeRNAsequencing/01_fastq/

nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr1.1.fastq.gz S906Nr1.2.fastq.gz S906Nr1.1.P.fastq.gz S906Nr1.1.U.fastq.gz S906Nr1.2.P.fastq.gz S906Nr1.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo1.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr2.1.fastq.gz S906Nr2.2.fastq.gz S906Nr2.1.P.fastq.gz S906Nr2.1.U.fastq.gz S906Nr2.2.P.fastq.gz S906Nr2.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo2.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr3.1.fastq.gz S906Nr3.2.fastq.gz S906Nr3.1.P.fastq.gz S906Nr3.1.U.fastq.gz S906Nr3.2.P.fastq.gz S906Nr3.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo3.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr4.1.fastq.gz S906Nr4.2.fastq.gz S906Nr4.1.P.fastq.gz S906Nr4.1.U.fastq.gz S906Nr4.2.P.fastq.gz S906Nr4.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo4.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr5.1.fastq.gz S906Nr5.2.fastq.gz S906Nr5.1.P.fastq.gz S906Nr5.1.U.fastq.gz S906Nr5.2.P.fastq.gz S906Nr5.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo5.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr6.1.fastq.gz S906Nr6.2.fastq.gz S906Nr6.1.P.fastq.gz S906Nr6.1.U.fastq.gz S906Nr6.2.P.fastq.gz S906Nr6.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo6.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr7.1.fastq.gz S906Nr7.2.fastq.gz S906Nr7.1.P.fastq.gz S906Nr7.1.U.fastq.gz S906Nr7.2.P.fastq.gz S906Nr7.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo7.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr8.1.fastq.gz S906Nr8.2.fastq.gz S906Nr8.1.P.fastq.gz S906Nr8.1.U.fastq.gz S906Nr8.2.P.fastq.gz S906Nr8.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo8.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr9.1.fastq.gz S906Nr9.2.fastq.gz S906Nr9.1.P.fastq.gz S906Nr9.1.U.fastq.gz S906Nr9.2.P.fastq.gz S906Nr9.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36       > Trimo9.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr10.1.fastq.gz S906Nr10.2.fastq.gz S906Nr10.1.P.fastq.gz S906Nr10.1.U.fastq.gz S906Nr10.2.P.fastq.gz S906Nr10.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > Trimo10.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr11.1.fastq.gz S906Nr11.2.fastq.gz S906Nr11.1.P.fastq.gz S906Nr11.1.U.fastq.gz S906Nr11.2.P.fastq.gz S906Nr11.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > Trimo11.log &
nohup java -Xmx16g -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar PE S906Nr12.1.fastq.gz S906Nr12.2.fastq.gz S906Nr12.1.P.fastq.gz S906Nr12.1.U.fastq.gz S906Nr12.2.P.fastq.gz S906Nr12.2.U.fastq.gz ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > Trimo12.log &

nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr1.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 1 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr1.1.P.fastq.gz S906Nr1.2.P.fastq.gz > tophat1.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr2.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 2 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr2.1.P.fastq.gz S906Nr2.2.P.fastq.gz > tophat2.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr3.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 3 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr3.1.P.fastq.gz S906Nr3.2.P.fastq.gz > tophat3.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr4.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 4 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr4.1.P.fastq.gz S906Nr4.2.P.fastq.gz > tophat4.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr5.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 5 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr5.1.P.fastq.gz S906Nr5.2.P.fastq.gz > tophat5.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr6.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 6 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr6.1.P.fastq.gz S906Nr6.2.P.fastq.gz > tophat6.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr7.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 7 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr7.1.P.fastq.gz S906Nr7.2.P.fastq.gz > tophat7.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr8.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 8 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr8.1.P.fastq.gz S906Nr8.2.P.fastq.gz > tophat8.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr9.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 9 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr9.1.P.fastq.gz S906Nr9.2.P.fastq.gz > tophat9.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr10.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 10 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr10.1.P.fastq.gz S906Nr10.2.P.fastq.gz > tophat10.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr11.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 11 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr11.1.P.fastq.gz S906Nr11.2.P.fastq.gz > tophat11.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr12.tophat2.aln.bam --transcriptome-index ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.95.genes/Sus_scrofa.Sscrofa11.1.95 --rg-id 906 --rg-sample 12 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/SScrofa/Sus_scrofa.Sscrofa11.1.dna.toplevel.bt2idx S906Nr12.1.P.fastq.gz S906Nr12.2.P.fastq.gz > tophat12.log &

# Just dump the unmapped reads
samtools bam2fq unmapped.bam > unmapped_reads.fastq
# Convert the reads from FASTQ to FASTA
sed -n '1~4s/^@/>/p;2~4p' unmapped_reads.fastq > unmapped_reads.fasta 
# Filter out the low complexity reads using the DUST algorithm
perl ~/Downloads/prinseq-lite-0.20.4/prinseq-lite.pl -fasta unmapped_reads.fasta -lc_method dust -lc_threshold 50
# Blast the Good unmapped reads to E Faecium.
blastn -task blastn -query unmapped_reads_prinseq_good_bBbZ.fasta -db ef.db -perc_identity 99 -outfmt 6 -evalue 0.1 -num_alignments 5 -out locations.txt


# When we want to dump the R1 and R2, for alignment against Efaecium, sort the bam file first
samtools sort -n unmapped.bam -o unmapped_sorted.bam
# Extract reads to two fastq files (ignore the reads that do not have a mate
/home/danny/Github/bedtools2/bin/bedtools bamtofastq -i unmapped_sorted.bam -fq unmapped_r1.fastq -fq2 unmapped_r2.fastq

### Prep the Efaecium genome for alignment
~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.fna.gz
~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.gff.gz

cd ~/genomes/EFaecium/
gunzip --keep GCF_000174395.2_ASM17439v2_genomic.fna.gz
gunzip --keep GCF_000174395.2_ASM17439v2_genomic.gff.gz

java -Xmx16g -jar /home/danny/picard-tools-1.119/CreateSequenceDictionary.jar R=~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.fna O=~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.dict
/home/danny/bowtie2-2.3.2/bowtie2-build -f ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.fna ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx
/home/danny/Github/samtools/samtools faidx ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.fna

ln -s GCF_000174395.2_ASM17439v2_genomic.fna GCF_000174395.2_ASM17439v2_genomic.bt2idx.fa

/home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -G GCF_000174395.2_ASM17439v2_genomic.gff --transcriptome-index=GCF_000174395.2_ASM17439v2_genomic.genes ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx

/home/danny/Github/STAR/source/STAR --runMode=genomeGenerate --genomeSAindexNbases 4 --genomeDir=/home/danny/genomes/EFaecium/STAR  --genomeFastaFiles /home/danny/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.fna --sjdbGTFfile /home/danny/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.gff --sjdbGTFtagExonParentTranscript Parent
/home/danny/Github/STAR/source/STAR --runMode alignReads --genomeDir=/home/danny/genomes/EFaecium/STAR --readFilesIn unmapped_reads.fastq --outSAMtype BAM SortedByCoordinate Unsorted


nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr1.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 1 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr2.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 2 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr3.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 3 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr4.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 4 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr5.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 5 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr6.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 6 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr7.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 7 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr8.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 8 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr9.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 9 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr10.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 10 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr11.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 11 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &
nohup /home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2 -o S906Nr12.EF.tophat2.aln.bam --transcriptome-index ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.genes/GCF_000174395.2_ASM17439v2_genomic --rg-id 906 --rg-sample 12 --rg-library RNA-seq --rg-platform Illumina --num-threads 8 ~/genomes/EFaecium/GCF_000174395.2_ASM17439v2_genomic.bt2idx unmapped_reads.fastq > tophatEF.log &

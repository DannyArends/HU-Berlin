# get the gffread utility
wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.3.Linux_x86_64.tar.gz

# Create a fasta file with all transcripts based on the reference genome, and the gene models in gff3
/home/danny/Github/gffread/gffread-0.12.3.Linux_x86_64/gffread -w alltranscripts.fa -g /home/danny/NAS/Cattle/Reference_Genomes/ARS-UCD1.2_Btau5.0.1Y/ARS-UCD1.2_Btau5.0.1Y.fa /home/danny/NAS/Cattle/Reference_Genomes/ARS-UCD1.2_Btau5.0.1Y/genes/Bos_taurus.ARS-UCD1.2.100.gff3

# run miranda using the bta-miR-200b
/home/danny/Downloads/miRanda-3.3a/src/miranda -strict miRNA.fa alltranscripts.fa -out miRanda.out

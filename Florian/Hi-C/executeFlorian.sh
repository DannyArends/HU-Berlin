# Digest the genome
~/Github/HU-Berlin/Florian/Hi-C/digestion/bin/DigestGenome > /home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68-digested.bed

# Iterative mapping BFMI (subset1)
nohup ~/Github/HU-Berlin/Florian/Hi-C/iterativeMapping/bin/IterativeMapping \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa \
/halde/Florian/raw_data/FKR02_1bfmi/subset/FKR02_1bfmi_R1_001.subset3.fastq.gz \
/halde/Florian/raw_data/FKR02_1bfmi/FKR02_1bfmi_R1.subset3.alignment > FKR02_1bfmi_R1.subset3.log 2>&1  &

nohup ~/Github/HU-Berlin/Florian/Hi-C/iterativeMapping/bin/IterativeMapping \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa \
/halde/Florian/raw_data/FKR02_1bfmi/subset/FKR02_1bfmi_R2_001.subset3.fastq.gz \
/halde/Florian/raw_data/FKR02_1bfmi/FKR02_1bfmi_R2.subset3.alignment > FKR02_1bfmi_R2.subset3.log 2>&1  &

# Assign fragments from both alignments (subset1)
~/Github/HU-Berlin/Florian/Hi-C/fragmentAssignment/bin/FragmentAssignment \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68-digested.bed \
/halde/Florian/raw_data/FKR02_1bfmi/FKR02_1bfmi_R1.subset1.alignment \
/halde/Florian/raw_data/FKR02_1bfmi/FKR02_1bfmi_R2.subset1.alignment \
/halde/Florian/raw_data/FKR02_1bfmi/FKR02_1bfmi.subset1.merged.alignments.txt

# Genome binning (subset1)
~/Github/HU-Berlin/Florian/Hi-C/genomeBinning/bin/GenomeBinning \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68-digested.bed \
/halde/Florian/raw_data/FKR02_1bfmi/FKR02_1bfmi.subset1.merged.alignments.txt \
/halde/Florian/raw_data/FKR02_1bfmi/FKR02_1bfmi.subset1.binned.alignments.txt



nohup ~/Github/HU-Berlin/Florian/Hi-C/iterativeMapping/bin/IterativeMapping \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa \
/halde/Florian/raw_data/FKR02_2b6n/subset/FKR02_2b6n_R1_001.subset3.fastq.gz \
/halde/Florian/raw_data/FKR02_2b6n/FKR02_2b6n_R1.subset3.alignment > FKR02_2b6n_R1.subset3.log 2>&1  &

nohup ~/Github/HU-Berlin/Florian/Hi-C/iterativeMapping/bin/IterativeMapping \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68.fa \
/halde/Florian/raw_data/FKR02_2b6n/subset/FKR02_2b6n_R2_001.subset3.fastq.gz \
/halde/Florian/raw_data/FKR02_2b6n/FKR02_2b6n_R2.subset3.alignment > FKR02_2b6n_R2.subset3.log 2>&1  &

# Assign fragments from both alignments
~/Github/HU-Berlin/Florian/Hi-C/fragmentAssignment/bin/FragmentAssignment \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68-digested.bed \
/halde/Florian/raw_data/FKR02_2b6n/FKR02_2b6n_R1.subset1.alignment \
/halde/Florian/raw_data/FKR02_2b6n/FKR02_2b6n_R2.subset1.alignment \
/halde/Florian/raw_data/FKR02_2b6n/FKR02_2b6n.subset1.merged.alignments.txt

# Genome binning
~/Github/HU-Berlin/Florian/Hi-C/genomeBinning/bin/GenomeBinning \
/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/GRCm38_68-digested.bed \
/halde/Florian/raw_data/FKR02_2b6n/FKR02_2b6n.subset1.merged.alignments.txt \
/halde/Florian/raw_data/FKR02_2b6n/FKR02_2b6n.subset1.binned.alignments.txt

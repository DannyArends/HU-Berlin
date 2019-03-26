
### Make a reference directory
mkdir /home/danny/References/Mouse/GRCm38_95
cd /home/danny/References/Mouse/GRCm38_95/

### Get the current reference
wget ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz

### Download trimmomatic and build from source
cd /home/danny/Github/trimmomatic/
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.38.zip
unzip Trimmomatic-Src-0.38.zip
cd trimmomatic-0.38/
ant

### Location of trimmomatic
#/home/danny/Github/trimmomatic/trimmomatic-0.38/dist/jar/trimmomatic-0.38.jar

### Download GATK
cd /home/danny/Github
wget https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
unzip gatk-4.1.0.0.zip

### Install picard tools
mkdir /home/danny/Github/picard-2.19.0
cd /home/danny/Github/picard-2.19.0
wget https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar

### BWA
cd /home/danny/Github/
git clone https://github.com/lh3/bwa.git
make -j 4

### Install samtools and htslib
cd /home/danny/Github/
git clone https://github.com/samtools/samtools.git
git clone https://github.com/samtools/htslib.git
cd samtools
make -j 4

### Create required index
cd /home/danny/References/Mouse/GRCm38_95
/home/danny/Github/bwa/bwa index Mus_musculus.GRCm38.dna.toplevel.fa.gz
/home/danny/Github/samtools/samtools faidx Mus_musculus.GRCm38.dna.toplevel.fa
java -jar /home/danny/Github/picard-2.19.0/picard.jar CreateSequenceDictionary REFERENCE=Mus_musculus.GRCm38.dna.toplevel.fa OUTPUT=Mus_musculus.GRCm38.dna.toplevel.dict 

### Make a folder to run the analysis
mkdir /halde/BFMI_Alignment_Mar19
cd /halde/BFMI_Alignment_Mar19

### Copy the pipeline script to /halde/BFMI_Alignment_Mar19
Rscript pipeline_Mar19.R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ext_L7258-1_BFMI860-S12_S1_R /halde/BFMI_Alignment_Mar19 BFMI860-12
Rscript pipeline_Mar19.R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ext_L7258-2_BFMI860-S12_S2_R /halde/BFMI_Alignment_Mar19 BFMI860-12
Rscript pipeline_Mar19.R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ext_L7258-1_BFMI860-S12_S9_R /halde/BFMI_Alignment_Mar19 BFMI860-12

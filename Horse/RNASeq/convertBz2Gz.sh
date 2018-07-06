bzip2 -dc 11N_R1.fastq.bz2 | gzip > 11N_R1.fastq.gz
bzip2 -dc 11V_R1.fastq.bz2 | gzip > 11V_R1.fastq.gz

bzip2 -dc 14N_R1.fastq.bz2 | gzip > 14N_R1.fastq.gz
bzip2 -dc 14V_R1.fastq.bz2 | gzip > 14V_R1.fastq.gz

bzip2 -dc 15N_R1.fastq.bz2 | gzip > 15N_R1.fastq.gz
bzip2 -dc 15V_R1.fastq.bz2 | gzip > 15V_R1.fastq.gz

bzip2 -dc 16N_R1.fastq.bz2 | gzip > 16N_R1.fastq.gz
bzip2 -dc 16V_R1.fastq.bz2 | gzip > 16V_R1.fastq.gz

bzip2 -dc 21N_R1.fastq.bz2 | gzip > 21N_R1.fastq.gz
bzip2 -dc 21V_R1.fastq.bz2 | gzip > 21V_R1.fastq.gz

bzip2 -dc 22N_R1.fastq.bz2 | gzip > 22N_R1.fastq.gz
bzip2 -dc 22V_R1.fastq.bz2 | gzip > 22V_R1.fastq.gz

bzip2 -dc 23N_R1.fastq.bz2 | gzip > 23N_R1.fastq.gz
bzip2 -dc 23V_R1.fastq.bz2 | gzip > 23V_R1.fastq.gz

bzip2 -dc 24N_R1.fastq.bz2 | gzip > 24N_R1.fastq.gz
bzip2 -dc 24V_R1.fastq.bz2 | gzip > 24V_R1.fastq.gz

nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/11V &> 11V.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/14V &> 14V.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/15V &> 15V.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/16V &> 16V.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/21V &> 21V.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/22V &> 22V.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/23V &> 23V.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/24V &> 24V.nohup.out &


nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/11N &> 11N.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/14N &> 14N.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/15N &> 15N.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/16N &> 16N.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/21N &> 21N.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/22N &> 22N.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/23N &> 23N.nohup.out &
nohup Rscript ~/RNAseqHorses/pipeline.R ~/NAS/Horse/RNA/RNAseq-Kabardian/24N &> 24N.nohup.out &








samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/AKR_J.bam 3:30000000-40000000 > AKR_J_Chr3.bam
samtools index AKR_J_Chr3.bam

samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZO_HlLtJ.bam 3:30000000-40000000 > NZO_HlLtJ_Chr3.bam
samtools index NZO_HlLtJ_Chr3.bam
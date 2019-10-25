samtools view -b /halde/MGP/ext_L7258-1_BFMI860-S12_S1_RP_trimmed.aligned.sorted.realigned.recalibrated.bam MT > BFMI_MT_1.bam
samtools view -b /halde/MGP/ext_L7258-2_BFMI860-S12_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam MT > BFMI_MT_2.bam
samtools merge BFMI_MT_merged.bam BFMI_MT_1.bam BFMI_MT_2.bam
samtools index BFMI_MT_merged.bam

samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZW_LacJ.bam MT > NZW_LacJ_MT.bam
samtools index NZW_LacJ_MT.bam

samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129P2_OlaHsd.bam MT > 129P2_OlaHsd_MT.bam
samtools index 129P2_OlaHsd_MT.bam

samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZO_HlLtJ.bam MT > NZO_HlLtJ_MT.bam
samtools index NZO_HlLtJ_MT.bam

samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam MT > C57BL_6NJ_MT.bam
samtools index C57BL_6NJ_MT.bam

samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/A_J.bam MT > A_J_MT.bam
samtools index A_J_MT.bam


samtools view -b /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SPRET_EiJ.bam MT > SPRET_EiJ_MT.bam
samtools index SPRET_EiJ_MT.bam
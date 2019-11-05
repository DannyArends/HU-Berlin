~/Github/bcftools/bcftools mpileup -r $1 -I -q 20 -Ou -f /home/danny/References/Mouse/GRCm38_68/GRCm38_68.fa \
/halde/MGP/ext_L7258-1_BFMI860-S12_S1_RP_trimmed.aligned.sorted.realigned.recalibrated.bam \
/halde/MGP/ext_L7258-2_BFMI860-S12_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129P2_OlaHsd.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/129S1_SvImJ.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/A_J.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/AKR_J.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C3H_HeJ.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/C57BL_6NJ.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/CAST_EiJ.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/DBA_2J.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZO_HlLtJ.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/NZW_LacJ.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/PWK_PhJ.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SJL_J.bam \
/home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/REL-1604-BAM/SPRET_EiJ.bam \
| ~/Github/bcftools/bcftools call -mv -Ou | ~/Github/bcftools/bcftools filter -i 'MIN(DP)>20 && QUAL>200' -Oz > MGP_BFMI_SNPs_Chr$1.vcf.gz

# Look at the header
gzip -cd Joint_69_flagged.tab.vcf.gz | head -n 570

# First index the CC vcf file
/home/danny/htslib/tabix -p vcf 

# Extract the regions we have
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr1:124548738-184926264 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL1.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr1:135464050-157672910 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL2.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr3:85146443-147514132 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL3.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr4:58831312-88222686 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL4.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr6:76940963-114257130 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL5.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr7:62305639-81923457 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL6.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr7:62305639-97260868 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL7.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr7:19482853-119720011 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL8.vcf
/home/danny/Github/bcftools/bcftools view -s CC011,CC017 -r chr8:31799796-104944836 Joint_69_flagged.tab.vcf.gz > Texas_Pb_QTL9.vcf

# SNP calling

### In case they do not provide a reference or DBSNP file, also download these files
# Get the reference 
wget ftp://ftp.ensembl.org/pub/release-87/fasta/gallus_gallus/dna/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz
# unzip the reference
gunzip Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz

# Get the known SNPs
wget ftp.ensembl.org/pub/release-87/variation/vcf/gallus_gallus/Gallus_gallus.vcf.gz

### Using the GATK
# For information on parameters and such: https://software.broadinstitute.org/gatk/gatkdocs/
# Download the GATK, from https://software.broadinstitute.org/gatk/download/
# Login: Danny.Arends@gmail.com Password: gatklogin
# Current version 3.7
# download to ~/tools/GATK
# extract: tar jxf GenomeAnalysisTK-3.7.tar.bz2

# SNP calling into a vcf file
java -Xmx12g -jar ~/tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller \                                  # GATK, using the haplotype caller
     -R ~/references/gallus_gallus5.0/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \                   # Genomic reference sequence
     -I ~/species1.bam -I ~/species2.bam -I ~/speciesN.bam \                                              # Input BAM files
     --dbsnp ~/references/gallus_gallus5.0/Gallus_gallus.vcf.gz \                                         # DBsnp reference
     -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 \                               # Settings
     -o mySNPs.vcf                                                                                        # Output file

# Variant filtration add information in the 4th column of the VCF file
java -Xmx4g -jar  ~/tools/GATK/GenomeAnalysisTK.jar -T VariantFiltration \                                # GATK, using the VariantFiltration
      -R ~/references/gallus_gallus5.0/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \                  # Genomic reference sequence
      -V mySNPs.vcf \                                                                                     # Input VCF file (Can come from GATK, or samtools)
      -window 35 -cluster 3 \                                                                             # Some parameters
      -filterName StrandBias -filter "FS > 30.0" \                                                        # Filter for strand bias
      -filterName LowQual -filter "QUAL < 30.0 || DP < 4" \                                               # Filter for low quality
      -o myFilteredSNPs.vcf                                                                               # Output file


### Call all 64 strains at the same time using the samtools
# For information on parameters and such: http://www.htslib.org/doc/samtools-0.1.19.html
~/tools/samtools/samtools mpileup \                                                                       # Mpileup using samtools
      -gf ~/references/gallus_gallus5.0/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \                 # Genomic reference sequence
       -o population.bcf \                                                                                # Output BCF file
       ~/species1.bam ~/species2.bam ~/speciesN.bam                                                       # Input BAM files

# Variant filtering step, removes SNPs failing the filter
~/tools/samtools/bcftools call -c population.bcf | ~/tools/samtools/vcfutils.pl varFilter -D 100 - > population.vcf
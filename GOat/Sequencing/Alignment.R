#Create indices
samtools faidx GoatGenome_6_8_14.fa
java -jar /home/shijie/tools/picard-tools-1.119/CreateSequenceDictionary.jar R=GoatGenome_6_8_14.fa O=GoatGenome_6_8_14.dict
/home/danny/Github/bwa/bwa-0.7.15/bwa index GoatGenome_6_8_14.fa

#Alignment

/home/danny/Github/bwa/bwa-0.7.15/bwa 
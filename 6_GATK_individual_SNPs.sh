#FIRST RESTRICT ALIGNMENTS TO UNIQUELY MAPPED READS (MAPq 255)   

#Melpomene adults
#remove first slurm jobs from /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/
cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load samtools/1.4.1'; echo "samtools view -b -q 255 dedupl.sorted.bam > unique.dedupl.sorted.bam") | sbatch
 done

#Cydno adults
cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load samtools/1.4.1'; echo "samtools view -b -q 255 dedupl.sorted.bam > unique.dedupl.sorted.bam") | sbatch
 done

#CREATE INDEXES FOR ALIGNMENT FILES
#Melpomene adults
cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load samtools/1.4.1'; echo "samtools index -b unique.dedupl.sorted.bam") | sbatch
 done

#Cydno adults
cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load samtools/1.4.1'; echo "samtools index -b unique.dedupl.sorted.bam") | sbatch
 done


#Split'N'Trim and reassign mapping qualities (#https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq)
#Use function SplitNCigarReads: splits reads into exon segments and hard-clip any sequences overhanging in the intronic regions

#(before that) need index and dictionary of the genome
cd /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/
samtools faidx Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa
#
cd /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/
java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar CreateSequenceDictionary R= Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa O= Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.dict

#Now Split and Trim
#Melpomene adults
cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP
  mkdir ./$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/$i/unique.dedupl.sorted.bam -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS") | sbatch
 done

#Cydno Adults 
cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP
  mkdir ./$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/$i/unique.dedupl.sorted.bam -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS") | sbatch
 done

#VARIANT CALLING
#Melpomene Adults     
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/output.vcf") | sbatch
 done

#Cydno Adults
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/output.vcf") | sbatch
 done

#COMMENTS
#actually GATK filters out not uniquely mapped reads anyway (<MAPQ20)
#-dontUseSoftClippedBases = take into account the information about intron-exon split regions that is embedded in the BAM file
#-stand_call_conf 20.0 = minimum phred-scaled confidence threshold for calling variants is 20
#by default the number of alternative alleles possible is 6 > maxAltAlleles (6)

#VARIANT FILTERING
#Melpomene Adults   
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/output.vcf -window 35 -cluster 3 -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/filtered.output.vcf") | sbatch
 done

#Cydno Adults   
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/output.vcf -window 35 -cluster 3 -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/filtered.output.vcf") | sbatch
 done

#COMMENTS
#how to read VCF files #https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it
#filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3
#filtering based on Fisher Strand values (FS > 30.0) and Qual By Depth values (QD < 2.0).

#Keep only filtered SNPs in the variant calling files (vcf)
#Melpomene Adults 
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/$i.filtered.output.vcf --excludeFiltered -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/$i.passed.output.vcf") | sbatch
 done

#Cydno Adults 
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/$i.filtered.output.vcf --excludeFiltered -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/$i.passed.output.vcf") | sbatch
 done

#repeat for other samples..

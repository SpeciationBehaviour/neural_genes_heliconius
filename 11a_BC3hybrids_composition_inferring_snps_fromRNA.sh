#JOINTLY GENOTYPE “Pure SPECIES” >>> to build database fixed MP and CP variants
#for joint genotyping need to run variant calling in different mode:
#https://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode
#and the jointly genotype:
#https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php



########156h
#Variant Calling
#Melpomene
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/$i/output.g.vcf") | sbatch
 done
#Cydno
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/output.g.vcf") | sbatch
 done
########60h
#Melpomene
cd /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/$i/output.g.vcf") | sbatch
 done
#Cydno
cd /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/$i/output.g.vcf") | sbatch
 done
######ADULTS
#Melpomene
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/output.g.vcf") | sbatch
 done
#Cydno
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/$i/output.g.vcf") | sbatch
 done



#CombineGVCFs all stages:
#https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php

#MP
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T CombineGVCFs -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/5/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/6/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/14/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/17/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/18/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/24/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/150/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/184/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/MP/220/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/87/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/92/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/95/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/97/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/115/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/117/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/124/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/149/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/164/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/MP/208/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/45/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/47/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/53/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/70/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/71A/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/78/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/80/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/83/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/100/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/104/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/128/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/218/output.g.vcf -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP/MP_combined.g.vcf") | sbatch

#CP
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T CombineGVCFs -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/8/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/13/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/21/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/30/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/137/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/142/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/151/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/156/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/168/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/85/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/86/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/90/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/118/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/119/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/125/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/144/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/146/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/162/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/CP/200/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/50/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/51/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/57/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/58/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/67/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/68/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/81/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/82/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/84/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/98/output.g.vcf --variant /data/home/wolfproj/wolfproj-06/7_GATK/Adults/CP/99/output.g.vcf -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/CP/CP_newcombined.g.vcf") | sbatch 


#removed one cydno individual #individual 4 --variant /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/4/output.g.vcf  #to use as control (see later)



#Combined genotyping
#GENOTYPE all stages MP - GenotypeGVCFs 
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP/MP_combined.g.vcf -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP/MP_genotyped.g.vcf") | sbatch
#it doesn’t include invariant sites

#GENOTYPE all stages CP - GenotypeGVCFs
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/CP/CP_newcombined.g.vcf -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/CP/CP_newgenotyped.g.vcf") | sbatch

 
#if want to extract only SNPS (not indels)
#MP
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP/MP_genotyped.g.vcf -selectType SNP -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP/MP_geno_SNPsonly.g.vcf") | sbatch
#CP
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/CP/CP_genotyped.g.vcf -selectType SNP -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/CP/CP_geno_SNPsonly.g.vcf") | sbatch




#FILTERING
#individuals filters: -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0'
#explanation parameters #https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it


#Melpomene 
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP/MP_genotyped.g.vcf -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_filtered_variants.vcf") | sbatch

#leave only SNPs/indels that PASSed the filter = with excludeFiltered. 
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_filtered_variants.vcf --excludeFiltered  -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_passed_variants.vcf") | sbatch

#now “filter” GENOTYPE FIELDS #set to null ./. #where DP < 4
#vcf tools ####https://www.biostars.org/p/200527/
(echo '#!/bin/bash'; echo '#SBATCH -J vcftools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14 '; echo "vcftools --vcf /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_passed_variants.vcf --out /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_geno_variants.vcf --minDP 4 --recode --recode-INFO-all") | sbatch

#remove where all genotypes null
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_geno_variants.vcf.recode.vcf --excludeNonVariants -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_variants.vcf") | sbatch


#Cydno 
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/CP/CP_newgenotyped.g.vcf -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_filtered_variants.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_filtered_variants.vcf --excludeFiltered -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_passed_variants.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J vcftools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14 '; echo "vcftools --vcf /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_passed_variants.vcf --out /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_geno_variants.vcf --minDP 4 --recode --recode-INFO-all") | sbatch

#
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_geno_variants.vcf.recode.vcf --excludeNonVariants -o /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_variants.vcf") | sbatch




####################################################
#CREATE DATABASE SNPs/indels MP and SNPs CP
#bgzip and create index
module load bcftools/1.4.1
bgzip MP_variants.vcf
bcftools index MP_variants.vcf.gz
#
bgzip CP_variants.vcf
bcftools index CP_variants.vcf.gz


#create files with genotypes unique to MP and genotypes unique to CP
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_variants.vcf.gz /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_variants.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/") | sbatch

#output files, first file (0000.vcf) contains variants that unique for the first file (MP). The second output (0001.vcf) contains variants unique for the second input file (CP). The third file contains variants intersected for both first and second input files (0002.vcf/ and 0003.vcf)


bgzip 0000.vcf
tabix -p vcf 0000.vcf.gz #MP
bgzip 0001.vcf
tabix -p vcf 0001.vcf.gz #CP



######## change coordinates from scaffolds to chromosomes
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/vcfChromTransfer.py /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/no_header_coordinates.txt /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/
#flip around coordinates # from scaffold to chromosomes
awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $5, $6, $1, $2, $3, $7}' no_header_coordinates.txt > coordinates.txt
#change header
newChrom        newStart        newEnd  chrom   start   end     orientation

#
(echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0000.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0000.chr.vcf") | sbatch

bgzip 0000.chr.vcf
tabix -p vcf 0000.chr.vcf.gz

#
(echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0001.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0001.chr.vcf") | sbatch

bgzip 0001.chr.vcf
tabix -p vcf 0001.chr.vcf.gz


#change coordinates also for hybrids vcfs
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/vcfChromTransfer.py /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/
#
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  (echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.output.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf") | sbatch
 done

#do the same for /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/E27
#then remove /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/vcfChromTransfer.py
#bzip files
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  bgzip /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf
 done
#
cd /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  bgzip /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/E27/$i/$i.passed.chr.output.vcf
 done
#then tabix
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  tabix -p vcf /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf.gz
 done
#
cd /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  tabix -p vcf /data/home/wolfproj/wolfproj-06/7_GATK/60h_APF/E27/$i/$i.passed.chr.output.vcf.gz
 done




###############################################
#INTERSECT 156h hybrids vcf with CP and MP SNPs database   #156h 
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP

#to CP
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0001.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i") | sbatch
 done


#to MP
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP

cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0000.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i") | sbatch
 done


#bgzip files #0003.vcf are intersects
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i/0003.vcf") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i/0003.vcf") | sbatch
 done




###SNP density/coverage 
#FIRST DIVIDE genome in INTERVALS/-kb windows
#copied Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa in /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows
#module load samtools/1.4.1
#samtools faidx Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa
#cut -f 1,2 Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa.fai > genome.sizes
3module load bedtools/2.26.0
#bedtools makewindows -g /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/genome.sizes -w 100000 > /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/windows.bed


#chromosomal coordinates
#nano new.genome.sizes
chr1    17206585
chr2    9045316
chr3    10541528
chr4    9662098
chr5    9908586
chr6    14054175
chr7    14308859
chr8    9320449
chr9    8708747
chr10   17965481
chr11   11759272
chr12   16327298
chr13   18127314
chr14   9174305
chr15   10235750
chr16   10083215
chr17   14773299
chr18   16803890
chr19   16399344
chr20   14871695
chr21   13359691

bedtools makewindows -g /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.genome.sizes -w 100000 > /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed




####### Calculate SNPs density #100kb windows #filtered variants
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i/CP_coverage.txt") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i/MP_coverage.txt") | sbatch
 done


#change names before downloading
#CP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toCP/$i/
  mv CP_coverage.txt $i.CP_coverage.txt    
 done
#MP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_filtered/toMP/$i/
  mv MP_coverage.txt $i.MP_coverage.txt    
 done

#download files
#go to R script




############# CONTROL WITH F1
   
# change coordinates
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/vcfChromTransfer.py /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/
#
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  (echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.passed.output.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.passed.chr.output.vcf") | sbatch
 done
#then remove /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/vcfChromTransfer.py

#bzip files
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  bgzip /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.passed.chr.output.vcf
 done
#then tabix
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  tabix -p vcf /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.passed.chr.output.vcf.gz
 done


#F1
#to CP
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0001.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i") | sbatch
 done

#to MP
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0000.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i") | sbatch
 done


#bgzip files
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i/0003.vcf") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i/0003.vcf") | sbatch
 done


#Calculate SNPs density #100kb windows
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i/CP_coverage.txt") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i/MP_coverage.txt") | sbatch
 done


#change names before downloading
#CP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toCP/$i/
  mv CP_coverage.txt $i.CP_coverage.txt    
 done
#MP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/F1_control_filtered/toMP/$i/
  mv MP_coverage.txt $i.MP_coverage.txt    
 done

#download files





###############CONTROL WITH PURE CYDNO
# change coordinates
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/vcfChromTransfer.py /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/
#
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  (echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/$i.passed.output.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/$i.passed.chr.output.vcf") | sbatch
 done
#then remove /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/vcfChromTransfer.py

#bzip files
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  bgzip /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/$i.passed.chr.output.vcf
 done
#then tabix
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  tabix -p vcf /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/$i.passed.chr.output.vcf.gz
 done



#156h CP
#INTERSECT 
#to CP
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0001.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i") | sbatch
 done
#to MP
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/0000.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/CP/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i") | sbatch
 done 


#bgzip
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i/0003.vcf") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i/0003.vcf") | sbatch
 done


#Calculate SNPs density #100kb windows
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i/CP_coverage.txt") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i/MP_coverage.txt") | sbatch
 done


#change names before downloading
#CP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toCP/$i/
  mv CP_coverage.txt $i.CP_coverage.txt    
 done
#MP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156_purespecies_filtered_control/toMP/$i/
  mv MP_coverage.txt $i.MP_coverage.txt    
 done
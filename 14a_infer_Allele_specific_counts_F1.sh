### FILTER VARIANTS INFERRED FROM GENOME RESEQUENCING DATA

#first bgzip files
cd /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqPASSED.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqPASSED.vcf") | sbatch

#tabix files
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "tabix /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqPASSED.vcf.gz") | sbatch 
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "tabix /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqPASSED.vcf.gz") | sbatch 

#SELECT ONLY SNPs (no indels)
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqPASSED.vcf.gz -selectType SNP -o /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_onlySNPs_allelefreqPASSED.vcf.gz") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqPASSED.vcf.gz -selectType SNP -o /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_onlySNPs_allelefreqPASSED.vcf.gz") | sbatch


## MAKE files with genotype positions unique to Melpomene and genotype positions unique to Cydno
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqPASSED.vcf.gz /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_onlySNPs_allelefreqPASSED.vcf.gz -p /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/") | sbatch
#bgzip
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0000.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0001.vcf") | sbatch
#tabix
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0000.vcf.gz") | sbatch #MP
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0001.vcf.gz") | sbatch #CP
#index 
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0000.vcf.gz") | sbatch 
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0001.vcf.gz") | sbatch 


## FILTER F1 hybrids vcfs files

#Variant Filtering 
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/output.vcf -filterName FS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' --filterName AF -filter 'AF > 0.6' -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.SNPcluster.filtered.output.vcf") | sbatch
 done

#keep only passed snps #and biallelic
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.SNPcluster.filtered.output.vcf --excludeFiltered --restrictAllelesTo BIALLELIC -o /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.SNPclustered.passed.output.vcf") | sbatch
 done

#bgzip files
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.SNPclustered.passed.output.vcf") | sbatch
 done
#index
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.SNPclustered.passed.output.vcf.gz") | sbatch
 done


## GET INTERSECT variants shared between F1 hybrids and cydno and melpomene (homozygous opposite in the species)

#1) variants inferred from cydno reads
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0001.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.SNPclustered.passed.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i") | sbatch
 done
#0003.vcf are intersects

#first sort vcf #order scaffolds must be the same #as in reference genome
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J sort'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar SortVcf I=/data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i/0003.vcf O=/data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i/0003.sort.vcf SEQUENCE_DICTIONARY=/data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.dict") | sbatch
 done

#then reindex files
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14'; echo "vcftools index /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i/0003.sort.vcf") | sbatch
 done

# ASE READER ##count how many reads map to cydno/melpomene allele
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_bam_files/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J ASE'; echo '#SBATCH -n 2'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T ASEReadCounter -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_bam_files/$i/unique.split.bam -sites /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/$i/0003.sort.vcf -U ALLOW_N_CIGAR_READS -minDepth 4 --minBaseQuality 2 -o /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/$i/allele_counts.csv") | sbatch
 done

#rename files before dowloading
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/$i/
  mv allele_counts.csv $i.allele_counts.csv    
 done


#2) variants inferred from melpomene reads
cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_allele_spec_expression_lessfiltered/0000.vcf.gz  /data/home/wolfproj/wolfproj-06/7_GATK/Adults/F1/$i/$i.SNPclustered.passed.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i") | sbatch
 done
#0003.vcf are intersects

#sort files
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J sort'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar SortVcf I=/data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i/0003.vcf O=/data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i/0003.sort.vcf SEQUENCE_DICTIONARY=/data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.dict") | sbatch
 done

#then reindex files
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14'; echo "vcftools index /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i/0003.sort.vcf") | sbatch
 done

## ASE READER ##count how many reads map to cydno/melpomene allele
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_bam_files/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J ASE'; echo '#SBATCH -n 2'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T ASEReadCounter -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_bam_files/$i/unique.split.bam -sites /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_hetero_vcf_files_withSNPclusters/toMP/$i/0003.sort.vcf -U ALLOW_N_CIGAR_READS -minDepth 4 --minBaseQuality 2 -o /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/toMP/$i/allele_counts.csv") | sbatch
 done

#rename files before downloading
cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/toMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/13_allele_specific_expression/F1_allele_counts_withSNPclusters/toMP/$i/
  mv allele_counts.csv $i.allele_counts.csv    
 done

#download files and go to next R script

#In case need CHANGE COORDINATES #chromosomal into Hmel2.5 scaffolds
#Simon’s script #https://github.com/simonhmartin/genomics_general/blob/master/VCF_processing/vcfChromTransfer.py
#cd /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files

#download Hmel2.5.chromosomes.agp. from Hmel2.5
#add chr in front of numbers (chromosomes) 
#remove scaffolds not placed

awk 'BEGIN {FS="\t"; OFS="\t"} {print $6, $7, $8, $1, $2, $3, $9}' Hmel2.5.chromosomes.agp.txt > coordinates.txt
#add new header 
#newChrom	newStart	newEnd	chrom	start	end orientation
#new name file “no_header_coordinates.txt”

(echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/ros10.Hmel2.bwa.default.HC.DP8.chromo.vcf.gz -t /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/no_header_coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/ros10.newcoord.Hmel2.5.bwa.default.HC.DP8.chromo.vcf") | sbatch

(echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/chi10.Hmel2.bwa.default.HC.DP8.chromo.vcf.gz -t /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/no_header_coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/chi10.newcoord.Hmel2.5.bwa.default.HC.DP8.chromo.vcf") | sbatch




######## Filtering genotype calls melpomene and cydno
#Rosina 
#“filter” GENOTYPE FIELDS #set to null ./. #where DP < 10, > 100  #GQ < 30
(echo '#!/bin/bash'; echo '#SBATCH -J vcftools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14 '; echo "vcftools --vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/ros10.Hmel2.bwa.default.HC.DP8.chromo.vcf.gz --out /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_filtered.vcf --minGQ 30 --minDP 10 --maxDP 100 --recode --recode-INFO-all") | sbatch

#remove where all genotypes null/ keep only variants sites
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_filtered.vcf.recode.vcf --excludeNonVariants -o /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_passed.vcf") | sbatch

#awk lines with PASS
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print }' ros10_geno_passed.vcf > ros10_geno_ONLYpassed.vcf


#Chioneus
(echo '#!/bin/bash'; echo '#SBATCH -J vcftools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14 '; echo "vcftools --vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/chi10.Hmel2.bwa.default.HC.DP8.chromo.vcf.gz --out /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_filtered.vcf --minGQ 30 --minDP 10 --maxDP 100 --recode --recode-INFO-all") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_filtered.vcf.recode.vcf --excludeNonVariants -o /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_passed.vcf") | sbatch
#
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print }' chi10_geno_passed.vcf > chi10_geno_ONLYpassed.vcf



cp /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_ONLYpassed.vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_ONLYpassed.vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/


#CREATE DATABASE SNPs/indels MP and SNPs CP
#bgzip #long..
cd /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/ros10_geno_ONLYpassed.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/chi10_geno_ONLYpassed.vcf") | sbatch

#and create index
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/ros10_geno_ONLYpassed.vcf.gz") | sbatch 
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/chi10_geno_ONLYpassed.vcf.gz") | sbatch 


#create files with genotypes unique to MP and genotypes unique to CP
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/ros10_geno_ONLYpassed.vcf.gz /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/chi10_geno_ONLYpassed.vcf.gz -p /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/") | sbatch


#
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.vcf.gz") | sbatch #MP
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.vcf.gz") | sbatch #CP



##### RETRANSFORM IN CHROMOSOMAL COORDINATES
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/vcfChromTransfer.py /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
#
python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf
python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf

(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf.gz") | sbatch #MP
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf.gz") | sbatch #CP



###### INTERSECT BC3 hybrids vcfs with CP and MP vcfs 
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP
#to CP
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i") | sbatch
 done


#to MP
#create directory /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i") | sbatch
 done



#bgzip files #0003.vcf are intersects
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/0003.vcf") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/0003.vcf") | sbatch
 done



#Calculate SNPs density #100kb windows #genomic variants
#to CP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/CP_coverage.txt") | sbatch
 done
#to MP
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/MP_coverage.txt") | sbatch
 done


#change names before downloading
#CP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/
  mv CP_coverage.txt $i.CP_coverage.txt    
 done
#MP coverage
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/
  mv MP_coverage.txt $i.MP_coverage.txt    
 done




#For CHANGING COORDINATES #Hmel2.5 chromosomal into scaffolds
#Simon’s script at: https://github.com/simonhmartin/genomics_general/blob/master/VCF_processing/vcfChromTransfer.py
#download Hmel2.5.chromosomes.agp. from Hmel2.5
#add chr in front of numbers (chromosomes) 
#remove scaffolds not placed
awk 'BEGIN {FS="\t"; OFS="\t"} {print $6, $7, $8, $1, $2, $3, $9}' Hmel2.5.chromosomes.agp.txt > coordinates.txt
#add new header: newChrom	newStart	newEnd	chrom	start	end orientation
#change coordinates
(echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/ros10.Hmel2.bwa.default.HC.DP8.chromo.vcf.gz -t /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/no_header_coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/ros10.newcoord.Hmel2.5.bwa.default.HC.DP8.chromo.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J coordinates '; echo '#SBATCH -n 1'; echo "python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/chi10.Hmel2.bwa.default.HC.DP8.chromo.vcf.gz -t /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/no_header_coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/chi10.newcoord.Hmel2.5.bwa.default.HC.DP8.chromo.vcf") | sbatch

#Filtering genotype calls in melpomene rosina and cydno chioneus (inferred from genome resequencing samples)
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

#copy files in another directory and work there
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_ONLYpassed.vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_ONLYpassed.vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/


#CREATE file with SNPs/indels fixed in melpomene and in cydno
#first bgzip files
cd /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/ros10_geno_ONLYpassed.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/chi10_geno_ONLYpassed.vcf") | sbatch
#and create index
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/ros10_geno_ONLYpassed.vcf.gz") | sbatch 
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/chi10_geno_ONLYpassed.vcf.gz") | sbatch 

#create files with genotypes unique to Melpomene and genotypes unique to Cydno
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/ros10_geno_ONLYpassed.vcf.gz /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/chi10_geno_ONLYpassed.vcf.gz -p /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/") | sbatch

#bgzip
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.vcf") | sbatch
#index
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.vcf.gz") | sbatch #MP
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.vcf.gz") | sbatch #CP

##### RETRANSFORM IN CHROMOSOMAL COORDINATES
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/vcfChromTransfer.py /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf
python vcfChromTransfer.py -v /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.vcf.gz -t /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/coordinates.txt >> /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf
#bgzip
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf") | sbatch
#index
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf.gz") | sbatch #MP
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf.gz") | sbatch #CP

###### GET INTERSECT variants between BC3 hybrids vcfs with Cydno and Melpomene vcfs 
#with cydno
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0001.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i") | sbatch
 done
#with melpomene
cd /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/0000.chr.vcf.gz /data/home/wolfproj/wolfproj-06/7_GATK/156h_APF/E27/$i/$i.passed.chr.output.vcf.gz -p /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i") | sbatch
 done

#0003.vcf are intersects
#bgzip files 
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/0003.vcf") | sbatch
 done
#
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bgzip /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/0003.vcf") | sbatch
 done

#Calculate number of SNPs density #in 100kb windows #that BC3 hyrbids share with:
#Cydno
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/CP_coverage.txt") | sbatch
 done
#Melpomene
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J bedtools'; echo '#SBATCH -n 1'; echo 'module load bedtools/2.26.0'; echo "bedtools coverage -a /data/home/wolfproj/wolfproj-06/10_introgression_line/genome_windows/new.windows.bed -b /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/0003.vcf.gz -counts > /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/MP_coverage.txt") | sbatch
 done

#change names before downloading
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toCP/$i/
  mv CP_coverage.txt $i.CP_coverage.txt    
 done
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/10_introgression_line/RESULTS/156h_genomic_data/toMP/$i/
  mv MP_coverage.txt $i.MP_coverage.txt    
 done

###########FILTERING SNPS FROM VCF


#1)

#Rosina 
#“filter” GENOTYPE FIELDS #set to null ./. #where DP < 10, > 100  #GQ < 30
(echo '#!/bin/bash'; echo '#SBATCH -J vcftools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14 '; echo "vcftools --vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/ros10.newcoord.Hmel2.5.bwa.default.HC.DP8.chromo.vcf --out /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_filtered.vcf --minGQ 30 --minDP 10 --maxDP 100 --recode --recode-INFO-all") | sbatch

#remove where all genotypes null/ keep only variants sites
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_filtered.vcf.recode.vcf --excludeNonVariants -o /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/ros10_geno_passed.vcf") | sbatch

#awk lines with PASS
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print }' ros10_geno_passed.vcf > ros10_geno_ONLYpassed.vcf

#Chioneus
(echo '#!/bin/bash'; echo '#SBATCH -J vcftools'; echo '#SBATCH -n 1'; echo 'module load vcftools/0.1.14 '; echo "vcftools --vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/vcf_files/chi10.newcoord.Hmel2.5.bwa.default.HC.DP8.chromo.vcf --out /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_filtered.vcf --minGQ 30 --minDP 10 --maxDP 100 --recode --recode-INFO-all") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_filtered.vcf.recode.vcf --excludeNonVariants -o /data/home/wolfproj/wolfproj-06/12_genomic_data/passed_snps/chi10_geno_passed.vcf") | sbatch
#
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print }' chi10_geno_passed.vcf > chi10_geno_ONLYpassed.vcf



#cp /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/chi10_geno_ONLYpassed.vcf.gz /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/
#cp /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/ros10_geno_ONLYpassed.vcf.gz /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/
#(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "gunzip /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_ONLYpassed.vcf.gz") | sbatch
#(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "gunzip /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_ONLYpassed.vcf.gz") | sbatch


#2)
 
#ADD AF 
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools +fill-tags /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_ONLYpassed.vcf -- -t AF") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 1'; echo 'module load bcftools/1.4.1'; echo "bcftools +fill-tags /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_ONLYpassed.vcf -- -t AF") | sbatch
#the slurm.* file contain the results # chi10_geno_withAF.vcf

#SELECT FOR AF > 0.9 
#remove biallelic (AF doesn’t work for multiple alleles)
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_withAF.vcf -filterName AF -filter 'AF < 0.9' -o /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqfiltered.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_withAF.vcf -filterName AF -filter 'AF < 0.9' -o /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqfiltered.vcf") | sbatch


########keep PASSED
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqfiltered.vcf --excludeFiltered --restrictAllelesTo BIALLELIC -o /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqPASSED.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqfiltered.vcf --excludeFiltered --restrictAllelesTo BIALLELIC -o /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqPASSED.vcf") | sbatch

 
#3)

#now need to filter for how many individuals have this variant in homozygous state #note that AF refers to frequency of alleles in samples that actually have the alternative allele 
#downloaded vcffilterjdk
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/FILTRATION_TOOL/jvarkit/dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.isHomVar()).count()> 4;' /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqPASSED.vcf > /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_4homozygous.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/FILTRATION_TOOL/jvarkit/dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.isHomVar()).count()> 4;' /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqPASSED.vcf > /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_4homozygous.vcf") | sbatch



cp /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_4homozygous.vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/
cp /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_4homozygous.vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/


#4)

#CREATE DATABASE SNPs/indels MP and SNPs CP
#bgzip
cd /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line/
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/ros10_4homozygous.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/chi10_4homozygous.vcf") | sbatch

#and create index
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/ros10_4homozygous.vcf.gz") | sbatch 
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools '; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools index /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/chi10_4homozygous.vcf.gz") | sbatch 


#create files with genotypes unique to MP and genotypes unique to CP
(echo '#!/bin/bash'; echo '#SBATCH -J bcftools'; echo '#SBATCH -n 5'; echo 'module load bcftools/1.4.1'; echo "bcftools isec /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/ros10_4homozygous.vcf.gz /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/chi10_4homozygous.vcf.gz -p /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/") | sbatch


#
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/0000.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J bgzip '; echo '#SBATCH -n 5'; echo "bgzip /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/0001.vcf") | sbatch
#
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/0000.vcf.gz") | sbatch #MP
(echo '#!/bin/bash'; echo '#SBATCH -J tabix '; echo '#SBATCH -n 5'; echo "tabix -p vcf /data/home/wolfproj/wolfproj-06/12_genomic_data/for_introgression_line_fixed/0001.vcf.gz") | sbatch #CP
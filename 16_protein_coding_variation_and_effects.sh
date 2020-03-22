#Filter variants inferred from RNAseq data #already filtered for FS, QD, DP 

#filter for AF > 0.8
#melpomene variants
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/CP_variants.vcf -filterName AF -filter 'AF < 0.8' -o /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_80_variants.vcf") | sbatch
#cydno variants
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/10_introgression_line/SNPs_database/MP_CP_variants_hardfiltered/MP_variants.vcf -filterName AF -filter 'AF < 0.8' -o /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_80_variants.vcf") | sbatch

#keep PASSED variants
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_80_variants.vcf --excludeFiltered -o /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_80_biaSNP_variants.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_80_variants.vcf --excludeFiltered -o /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_80_biaSNP_variants.vcf") | sbatch

#keep variants found in at least 7 individuals 
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/FILTRATION_TOOL/jvarkit/dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.isHomVar()).count()> 7;' /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_80_biaSNP_variants.vcf > /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_AF80_7ind_variants.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/FILTRATION_TOOL/jvarkit/dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.isHomVar()).count()> 7;' /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_80_biaSNP_variants.vcf > /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_AF80_7ind_variants.vcf") | sbatch

#keep only opposite alleles
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_AF80_7ind_variants.vcf
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_AF80_7ind_variants.vcf

#create indexes
module load bcftools/1.4.1
bcftools index /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_AF80_7ind_variants.vcf.gz
bcftools index /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_AF80_7ind_variants.vcf.gz

#create files with genotypes unique to Melpomene and genotypes unique to Cyndo
bcftools isec /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/MP_AF80_7ind_variants.vcf.gz /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/CP_AF80_7ind_variants.vcf.gz -p /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/

#bzip files
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/0000.vcf
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/0001.vcf
#tabix files
tabix -p vcf /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/0000.vcf.gz #MP
tabix -p vcf /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/0001.vcf.gz #CP



#FILTER VARIANTS INFERRED FROM GENOME RESEQENCING DATA DATA (already filtered) for DP, GQ and AF>0.9

#keep variants that passed the filters
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/ros10_geno_allelefreqfiltered.vcf --excludeFiltered -o /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/ros10_geno_allelefreqPASSED.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/7_GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V /data/home/wolfproj/wolfproj-06/12_genomic_data/fixed_snps/chi10_geno_allelefreqfiltered.vcf --excludeFiltered -o /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/chi10_geno_allelefreqPASSED.vcf") | sbatch

#select those homozygous variants found in at least 8 individuals
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/FILTRATION_TOOL/jvarkit/dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.isHomVar()).count()> 8;' /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/ros10_geno_allelefreqPASSED.vcf > /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/ros10_8homozygous.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J GATK'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/FILTRATION_TOOL/jvarkit/dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.isHomVar()).count()> 8;' /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/chi10_geno_allelefreqPASSED.vcf >/data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/chi10_8homozygous.vcf") | sbatch

#keep only opposite alleles
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/ros10_8homozygous.vcf
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/chi10_8homozygous.vcf

#create index
module load bcftools/1.4.1
bcftools index /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/ros10_8homozygous.vcf.gz
bcftools index /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/chi10_8homozygous.vcf.gz

#create files with genotypes unique to MP and genotypes unique to CP
bcftools isec /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/ros10_8homozygous.vcf.gz /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/chi10_8homozygous.vcf.gz -p /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/
#
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/0000.vcf
bgzip /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/0001.vcf
#
tabix -p vcf /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/0000.vcf.gz #MP
tabix -p vcf /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/0001.vcf.gz #CP



### INTERSECT THE TWO DATASETS (VARIANTS INFERRED FROM RNA-seq DATA and from GENOMER-RESEQUENCING DATA
/data/home/wolfproj/wolfproj-06/15_protein_coding_changes/intersectvariants

#intersect melpomene substitutions
bcftools isec /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/0000.vcf.gz /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/0000.vcf.gz -p /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/intersectvariants/MP/
#intersext cydno substitutions
bcftools isec /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/RNAvariants/0001.vcf.gz /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/genomicvariants/0001.vcf.gz -p /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/intersectvariants/CP/



### ANNOTATE VARIANTS EFFECTS with SnpEff (http://snpeff.sourceforge.net/SnpEff_manual.html)

#FIRST HAVE TO BUILD DATABASE HELICONIUS GENOME  (http://snpeff.sourceforge.net/SnpEff_manual.html#databases) 
#http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_ENSEMBL_BFMPP_32_161.zip  
#add genome to configuration file:
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/snpEff.config
nano  snpEff.config
#add below 'Mouse / mm9.genome : Mouse / mm10.genome : Mouse', the following lines:
# Heliconius genome, version Hmel2.5
# Hmel2.5.genome : Heliconius
mkdir /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/data/Hmel2.5
cp /data/home/wolfproj/wolfproj-06/Genome_annotations/Hmel2.5.gff3 /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/data/Hmel2.5/ 
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/data/Hmel2.5/
mv Hmel2.5.gff3 genes.gff
cp /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/data/Hmel2.5/ 
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/data/Hmel2.5/
mv Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa sequences.fa
#build database
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff
java -jar snpEff.jar build -gff3 -v Hmel2.5 
cd /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/data/Hmel2.5


## NOW filter for SNPs that are found within protein-coding region
(echo '#!/bin/bash'; echo '#SBATCH -J SNPeff'; echo '#SBATCH -n 1'; echo "java -Xmx4g -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/snpEff.jar -canon -onlyProtein -no-downstream -no-intergenic -no-intron -no-upstream -no-utr Hmel2.5 /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/intersectvariants/MP/0003.vcf > /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/MP/MPfixed_ann.vcf") | sbatch
(echo '#!/bin/bash'; echo '#SBATCH -J SNPeff'; echo '#SBATCH -n 1'; echo "java -Xmx4g -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/snpEff.jar -canon -onlyProtein -no-downstream -no-intergenic -no-intron -no-upstream -no-utr Hmel2.5 /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/intersectvariants/CP/0003.vcf > /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/CP/CPfixed_ann.vcf") | sbatch

#COMMENTS
#canon : Only use canonical transcripts #Canonical transcripts are defined as the longest CDS of amongst the protein coding transcripts in a gene. If none of the transcripts in a gene is protein coding, then it is the longest cDNA.
#onlyProtein > use only protein coding transcripts #ANN field 


## SELECT ONLY MODERATE OR HIGH IMPACT VARIANTS
java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/SnpSift.jar filter "ANN[*].IMPACT has 'HIGH' | ANN[*].IMPACT has 'MODERATE'" /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/MP/MPfixed_ann.vcf > /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/MP/MPfixed_withimpact.vcf
java -jar /data/home/wolfproj/wolfproj-06/10_introgression_line/snpEff/SnpSift.jar filter "ANN[*].IMPACT has 'HIGH' | ANN[*].IMPACT has 'MODERATE'" /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/CP/CPfixed_ann.vcf > /data/home/wolfproj/wolfproj-06/15_protein_coding_changes/CP/CPfixed_withimpact.vcf








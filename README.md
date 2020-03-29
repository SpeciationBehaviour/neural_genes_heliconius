## Analysis details for “Visual mate preference evolution during butterfly speciation is linked to neural processing genes” ##

Matteo Rossi, Alexander E. Hausmann, Timothy J. Thurman,  Stephen H. Montgomery, Riccardo Papa, Chris D. Jiggins , W. Owen McMillan & Richard M. Merrill

This repository documents analyses for the manuscript "Visual mate preference evolution during butterfly speciation is linked to neural processing genes". These bespoke scripts, written for either Unix or R, are presented for transparency only and will require editing to be applied to other datasets. 

The repository includes:

### i) Raw data for behavioral analysis:
`ternary_final.csv`

### ii) A markdown file containing details for analyses of male courtship initiation data:
`0_analyses_behavioral_data.md`

### iii) Scripts for analyses of NGS data

1. Run a fast quality check on RNA-seq raw reads:
`1_script_fastqc.sh`

2. Trim adaptor and low quality bases with TrimGalore (wrapper script of Cutadapt):
`2_script_trimgalore.sh`

3. Map trimmed reads to the Heliconius melpomene (Hmel2.5) or H. cydno genomes with STAR:
`3_script_STAR.sh`

4. Count for each sample the number of RNA reads mapping to each gene model: 
`4_script_HTSeq.sh` 

5. Mark duplicate RNA reads (first step for variant calling in GATK, next script):
`5_script_Picard_markduplicates.sh`

6. Call variants with GATK from individual RNA samples:
`6_GATK_individual_SNPs.sh`

7. Estimate heterozygosity percentages on the Z chromosome to sex pupae:
`7_Sexing.R`

8. conduct differential gene expression analyses between pure species and hybrids with DESeq2 (gene counts estimated in 4_script_HTSeq):
`8_Differential_expression.R`

9. Use InterProScan and Panther to retrieve annotated functions associated with protein sequences of Hmel2.5: 
`9_InterProScan_Panther.sh` 

11. a) Infer BC3 hybrids genome composition, using variants inferred from RNA-seq data:
`11a_BC3hybrids_composition_inferring_snps_fromRNA.sh` 

b) Infer BC3 hybrids genome composition, using variants inferred from genome resequencing data:
`11b_BC3hybrids_composition_inferring_snps_from_genome_resequencing_data.sh`

12. Conduct differential gene expression analyses between BC3 hybrids segregating at the QTL for behavior on chromosome 18:
`12_Differential_expression_introgressionline.R`

13.  Estimate variants fixed in H. melpomene and H. cydno from genome resequencing data (for allele specific expression analyses):
`13_infer_fixed_snps_from_genomic_data.sh` 

14. a) Count RNA reads in F1 hybrids that map to either the cydno or the melpomene allele:
`14a_infer_Allele_specific_counts_F1.sh`

b) Conduct allele-specific expression analyses in hybrids + make Figure 3:
`14b_allele_specific_differential_expression.R`

15. a) Reannotate Hmel2.5 gene models using a transcript-based approach with Cufflinks:
`15a_reannotate_transcript_based_Cufflinks.sh`

b) Map trimmed reads to the H. melpomene genome and count reads mapping to gene models as annotated with Cufflinks in the previous script:
`15b_remapandcount_annotationCufflinks.sh`

c) Conduct differential gene expression analyses between pure species using the transcript-based annotated gene models: 
`15c_differential_expression_with_Cufflinks_annotation.sh`

16. Estimate variants fixed in H. melpomene and H. cydno from both RNAseq and genome resequencing data, and annotate variants effects with SNPeff:
`16_protein_coding_variation_and_effects.sh`

### iv) Scripts or markdown files to generate figures

`Figure_1_and_supporting_figure_1_preference_behaviour.md`

`Figure2_differential_expression_all_stages_combined.R`

`Figure4_admixture_graph.R`

`Supporting_figure_2_mapping_to_melpomene/cydno_genomes.R`

`Supporting_figure_5B_BC3hybrids_genome_composition.R`

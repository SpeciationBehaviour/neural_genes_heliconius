#Before mapping RNA-seq reads to reference genome with STAR, need to generate (reference) genome index
#Hmel2.5 reference and annotation downloaded at: http://butterflygenome.org/?q=node/4
#Hcydno reference and annotation can be found in the future as supplement from "K. Darragh, et al., A novel terpene synthase produces an anti-aphrodisiac pheromone in the butterfly Heliconius melpomene. bioRxiv, 779678 (2019)."
#in the annotations (when using STAR) tabs are not allowed in chromosomes/scaffold names and spaces are not recommended
#change header in  Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa from  “>Hmel201001o heliconius_melpomene_melpomene_hmel25_core_32_85_1 scaffold” to “>Hmel201001o” 
sed 's/ .*$//' Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa> newHeliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa 

#MELPOMENE genome index generation  
#cd /data/home/wolfproj/wolfproj-06/3_STAR/
(echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Hmel2_index/ --genomeFastaFiles /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --sjdbGTFfile /data/home/wolfproj/wolfproj-06/Genome_annotations/Hmel2.5.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99") | sbatch

#CYDNO genome index generation 
#(see scripts 11_ for correction of gene names)
(echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Cyd_index/ --genomeFastaFiles /data/home/wolfproj/wolfproj-06/Genome_assemblies/Cydno/cydRagout_correct.fa --sjdbGTFfile /data/home/wolfproj/wolfproj-06/Genome_annotations/Hcyd_corrected.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99") | sbatch

#NOTES
#—genomeFastaFiles specifies one or more FASTA files with the genome reference sequences #multiple reference sequences (henceforth called chromosomes) are allowed for each fasta file  #it is advised to add unplaced and unlocalized scaffolds 
#—sjdbGTFfile specifies the path to the file with annotated transcripts in the standard GTF format. STAR will extract splice junctions from this file #Chromosome names in the annotations GTF file have to match chromosome names in the FASTA genome sequence files. 
#for GFF3 formatted annotations you need to use --sjdbGTFtagExonParentTranscript Parent. For --sjdbGTFfile files STAR only processes lines which have --sjdbGTFfeatureExon  in the 3rd field (column) (so =exon by default). The exons are assigned to the transcripts using parent-child relationship defined by the --sjdbGTFtagExonParentTranscript (=transcript id by default) in GTF attribute  
#—sjdbOverhang specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina 2x100b paired-end reads, the ideal value is 100-1=99. ##In case of reads of varying length, the ideal value is max(ReadLength)-1. #set to 99 as we have 100bp paired-end reads


#MAPPING reads with STAR, in 2-pass mode

#1-pass
#Map reads from melpomene adults to the melpomene genome  (do for every file listed in a given directory)
#first create directory: /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoMP
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoMP
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoMP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Hmel2_index/ --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_2.fq.gz --readFilesCommand zcat --outFilterType BySJout --outSJfilterIntronMaxVsReadN 50000 100000 150000 --outSAMattributes NH HI AS nM XS") | sbatch
 done
 
##Map reads from cydno adults to the melpomene genome
#first create directory: /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/CPtoMP
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/CPtoMP
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/CPtoMP/$i/ 
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Hmel2_index/ --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_2.fq.gz --readFilesCommand zcat --outFilterType BySJout --outSJfilterIntronMaxVsReadN 50000 100000 150000 --outSAMattributes NH HI AS nM XS") | sbatch
 done

#last example #Map reads from melpomene adults to the cydno genome   
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoCP
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoCP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Cyd_index/ --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_2.fq.gz --readFilesCommand zcat --outFilterType BySJout --outSJfilterIntronMaxVsReadN 50000 100000 150000 --outSAMattributes NH HI AS nM XS") | sbatch
 done

#NOTES 
# —outSJfilterIntronMaxVsReadN 50000 100000 150000 #set last to 150000 because it’s (a bit more than) longest intron length detected in the Hmel2.5 annotation


#2nd-pass
#Mel to Mel
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Hmel2_index/ --sjdbFileChrStartEnd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoMP/*/*J.out.tab --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_2.fq.gz --readFilesCommand zcat --outSAMattributes NH HI AS nM XS") | sbatch
 done

#Cyd to Mel
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Hmel2_index/ --sjdbFileChrStartEnd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/CPtoMP/*/*J.out.tab --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_2.fq.gz --readFilesCommand zcat --outSAMattributes NH HI AS nM XS") | sbatch
 done

#Mel to Cyd
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCP
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCP/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/Cyd_index/ --sjdbFileChrStartEnd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoCP/*/*J.out.tab --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_2.fq.gz --readFilesCommand zcat --outSAMattributes NH HI AS nM XS") | sbatch
 done


#NOTES 
#the second pass does not increase the number of detected novel junctions, but allows to detect more spliced reads mapping to novel junctions. Run 1st pass of STAR mapping with the usual parameters, then collect the junctions detected in the first pass, and use them as "annotated" junctions for the 2nd pass mapping #if using annotation, discover new splice junctions then these junctions are merged with the GFF file #The filtering (and 2nd-pass gain) is only done for newly discovered junctions


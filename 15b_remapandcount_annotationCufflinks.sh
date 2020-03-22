#Map RNA-reads to the newly created transcript based annotations of the melpomene genome
#Example with Adult samples
#first generate genome index (with transcript based annotation) 
cd /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults
(echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults --genomeFastaFiles /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --sjdbGTFfile /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/combined_newpar13/combined_newpar13_stats/transcripts.gtf --sjdbOverhang 99") | sbatch
#removed --sjdbGTFtagExonParentTranscript Parent #only for gff format

#Use STAR in 2-pass mode

#1-pass
#Melpomene reads
#before create directory: /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoCuff
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoCuff
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoCuff/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_2.fq.gz --readFilesCommand zcat --outFilterType BySJout --outSJfilterIntronMaxVsReadN 50000 100000 150000 --outSAMattributes NH HI AS nM XS") | sbatch
 done

#cydno reads
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/CPtoCuff/
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/CPtoCuff/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_2.fq.gz --readFilesCommand zcat --outFilterType BySJout --outSJfilterIntronMaxVsReadN 50000 100000 150000 --outSAMattributes NH HI AS nM XS") | sbatch
 done

#F1 hybrids reads
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/F1toCuff/
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/F1toCuff/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/F1/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/F1/$i/*_2.fq.gz --readFilesCommand zcat --outFilterType BySJout --outSJfilterIntronMaxVsReadN 50000 100000 150000 --outSAMattributes NH HI AS nM XS") | sbatch
 done

#2-pass
#Melpomene
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCuff
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCuff/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults --sjdbFileChrStartEnd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/MPtoCuff/*/*J.out.tab --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/*_2.fq.gz --readFilesCommand zcat --outSAMattributes NH HI AS nM XS") | sbatch
 done
#Cydno
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoCuff
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoCuff/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults --sjdbFileChrStartEnd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/CPtoCuff/*/*J.out.tab --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/*_2.fq.gz --readFilesCommand zcat --outSAMattributes NH HI AS nM XS") | sbatch
 done
#F1 hybrids
cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/F1
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toCuff
  mkdir ./$i/
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toCuff/$i/
  (echo '#!/bin/bash'; echo '#SBATCH -J STAR'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load STAR/2.4.2a'; echo "STAR --runThreadN 1 --genomeDir /data/home/wolfproj/wolfproj-06/3_STAR/cufflink_index/Adults --sjdbFileChrStartEnd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/1stPass/F1toCuff/*/*J.out.tab --readFilesIn /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/F1/$i/*_1.fq.gz /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/F1/$i/*_2.fq.gz --readFilesCommand zcat --outSAMattributes NH HI AS nM XS") | sbatch
 done


#Before counting reads mapping to new gene models
#KEEP ONLY PROPERLT MAPPED READ PAIRS with Samtools

#Melpomene read alignments
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCuff/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCuff/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load samtools/1.4.1'; echo "samtools view -h -f 0x02 Aligned.out.sam > Aligned.out.proper.sam") | sbatch
 done
#Cydno read alignments
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoCuff/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoCuff/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load samtools/1.4.1'; echo "samtools view -h -f 0x02 Aligned.out.sam > Aligned.out.proper.sam") | sbatch
 done
#F1 hybrids read alignments
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toCuff/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toCuff/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load samtools/1.4.1'; echo "samtools view -h -f 0x02 Aligned.out.sam > Aligned.out.proper.sam") | sbatch
 done


#COUNT reads mapping to new gene models in Melpomene samples
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCuff/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/MPtoCuff/
  (echo '#!/bin/bash'; echo '#SBATCH -J HTSeq'; echo '#SBATCH -n 1';  echo '#SBATCH -t 23:59:59'; echo 'module load htseq/0.9.1'; echo "htseq-count -r pos -a 20 -m union --stranded=no -t exon --idattr gene_id /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoCuff/$i/Aligned.out.proper.sam /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/combined_newpar13/combined_newpar13_stats/transcripts.gtf > $i.propergenecounts.txt") | sbatch
 done
#make table for gene counts
cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/MPtoCuff/
#remove slurm.out
FILES=$(ls -t -v *.txt | tr '\n' ' ');
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > merged_htseq_MPtoCuff_genecounts.txt
#change header of merged_htseq_MPtoCuff_genecounts.txt in:
#gene_id	MP45_A_m	MP47_A_m	MP53_A_f  MP70_A_m	MP71_A_m	MP78_A_f 	MP80_A_f 	MP83_A_m	MP100_A_m	MP104_A_m	MP128_A_f	MP218_A_f
#remove these rows from merged_htseq_MPtoCuff_genecounts.txt	
__too_low_aQual	0	0	0	0	0	0	0	0	0	0	0	0
__ambiguous	2184175	1118184	1992730	1811619	1805012	1614617	1515048	1530842	2312038	2377535	2119335	2316873
__no_feature	482135	242688	442274	401624	357338	218519	240897	345653	606762	700699	544758	423891
__alignment_not_unique	861179	350659	854972	688237	662725	704423	563901	566903	878690	876712	783562	971602
__not_aligned	0	0	0	0	0	0	0	0	0	0	0	0


#COUNT reads mapping to new gene models in Cydno samples 
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoCuff/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/CPtoCuff/
  (echo '#!/bin/bash'; echo '#SBATCH -J HTSeq'; echo '#SBATCH -n 1';  echo '#SBATCH -t 23:59:59'; echo 'module load htseq/0.9.1'; echo "htseq-count -r pos -a 20 -m union --stranded=no -t exon --idattr gene_id /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoCuff/$i/Aligned.out.proper.sam /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/combined_newpar13/combined_newpar13_stats/transcripts.gtf > $i.propergenecounts.txt") | sbatch
 done
#make table
cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/CPtoCuff/
#remove slurm.out
FILES=$(ls -t -v *.txt | tr '\n' ' ');
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > merged_htseq_CPtoCuff_genecounts.txt
#change header: 
#gene_id CP50_A_f	CP51_A_f	CP57_A_m	CP58_A_f	CP67_A_f	CP68_A_f	CP81_A_f	CP82_A_m	CP84_A_f	CP98_A_m	CP99_A_m
#remove from file:	
__too_low_aQual	0	0	0	0	0	0	0	0	0	0	0
__ambiguous	2579829	2190832	1381838	2678079	1523859	968415	1382839	1367290	2008239	1932016	1826951
__no_feature	445018	409321	307129	547052	288161	193083	204193	409216	388777	399256	492378
__alignment_not_unique	1070548	905314	580331	1082689	545638	391973	537123	541481	783883	670116	728817


#COUNT reads mapping to new gene models in F1 hybrids
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toCuff/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/F1toCuff/
  (echo '#!/bin/bash'; echo '#SBATCH -J HTSeq'; echo '#SBATCH -n 1';  echo '#SBATCH -t 23:59:59'; echo 'module load htseq/0.9.1'; echo "htseq-count -r pos -a 20 -m union --stranded=no -t exon --idattr gene_id /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toCuff/$i/Aligned.out.proper.sam /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/combined_newpar13/combined_newpar13_stats/transcripts.gtf > $i.propergenecounts.txt") | sbatch
 done

#make table
cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/F1toCuff/
#remove slurm.out
FILES=$(ls -t -v *.txt | tr '\n' ' ');
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > merged_htseq_F1toCuff_genecounts.txt
#change header merged_htseq_F1toCuff_genecounts.txt:
#gene_id	F1_42_m	F1_48_m	F1_49_m	F1_56_f	F1_69_f
#remove these rows from merged_htseq_F1toCuff_genecounts.txt:	
__too_low_aQual	0	0	0	0	0
__ambiguous	1189578	1475189	1414186	1831805	1645884
__no_feature	232214	295220	357470	388982	331531
__alignment_not_unique	371847	525588	382782	556376	584750

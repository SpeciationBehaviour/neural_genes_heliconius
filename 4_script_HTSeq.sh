#Before counting RNA-reads that mapped to each gene model with HTseq, run Samtools to keep only "properly-mapped" pairs
#Example, keep only "properly-mapped" pairs for mapping files of every melpomene adult sample (mapped to the melpomene genome)
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load samtools/1.4.1'; echo "samtools view -h -f 0x02 Aligned.out.sam > Aligned.out.proper.sam") | sbatch
 done

#COMMENTS
#-h include the header
#-f 0x02 the read is mapped in a proper pair #in this way it filters out also #f 0x08 = if the mate is unmapped

#other example #cydno adults files
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load samtools/1.4.1'; echo "samtools view -h -f 0x02 Aligned.out.sam > Aligned.out.proper.sam") | sbatch
 done

#F1 hybrids files
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/F1toMP/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J samtools'; echo '#SBATCH -n 1'; echo 'module load samtools/1.4.1'; echo "samtools view -h -f 0x02 Aligned.out.sam > Aligned.out.proper.sam") | sbatch
 done


#Run HTseq on mapping files in ".sam" format
#example for files of every melpomene adult
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/MPtoMP/
  (echo '#!/bin/bash'; echo '#SBATCH -J HTSeq'; echo '#SBATCH -n 1'; echo 'module load htseq/0.9.1'; echo "htseq-count -r pos -a 20 -m union --stranded=no -t gene --idattr ID /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/$i/Aligned.out.proper.sam /data/home/wolfproj/wolfproj-06/Genome_annotations/Hmel2.5.gff3 > $i.propergenecounts.txt") | sbatch
 done

#NOTE
#-a = mapping quality filtering of the alignment (>20)

#Make table containg gene counts for every melpomene adult sample
cd /data/home/wolfproj/wolfproj-06/5_HTseq/Adults/MPtoMP/
#remove "*slurm.out" files
FILES=$(ls -t -v *.txt | tr '\n' ' ');
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > merged_htseq_MPtoMP_genecounts.txt

#COPY FILES on laptop and format it: 

#add header #for example for the file "merged_htseq_MPtoMP_genecounts.txt":
#gene_id	MP45_A_m	MP47_A_m	MP53_A_f  MP70_A_m	MP71_A_m	MP78_A_f 	MP80_A_f 	MP83_A_m	MP100_A_m	MP104_A_m	MP128_A_f	MP218_A_f

#then removed from file the following rows:	
#__no_feature	3020314	1655011	2780154	2376162	2362679	1887295	1837890	2147028	3393217	3576709	3091746	3010553
#__ambiguous	66906	24835	60686	67677	69527	89765	75985	48658	80402	79752	71152	94127 
#__too_low_aQual	0	0	0	0	0	0	0	0	0	0	0	0                                     
#__not_aligned	0	0	0	0	0	0	0	0	0	0	0	0
#__alignment_not_unique	863129	350828	856688	689523	664363	707609	565544	567350	882234	878921	785625	973952

#repeat the above procedure for all other species or hybrids files

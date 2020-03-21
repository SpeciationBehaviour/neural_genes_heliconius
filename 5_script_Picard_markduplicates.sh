#needs sorted BAM #transform SAM files with Picard
#Adults
#Melpomene (to MP)
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  (echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar SortSam I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/$i/Aligned.out.proper.sam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/$i/sorted.Aligned.out.bam SORT_ORDER=coordinate") | sbatch
 done
#then remove slurm jobs from directory
#Cydno (to MP)
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  (echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar SortSam I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/$i/Aligned.out.proper.sam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/$i/sorted.Aligned.out.bam SORT_ORDER=coordinate") | sbatch
 done
#then remove slurm jobs from directory


##Add read group information 
#instrument name:run id:flowcell id:lane

MP 45 HWI-D248:132:H8LHJADXX:2
MP 47 HISEQ:192:H9C0LADXX:1           
MP 53 HWI-D248:132:H8LHJADXX:2
MP 70 HWI-D248:132:H8LHJADXX:2
MP 71 HISEQ:192:H9C0LADXX:1           
MP 78 HWI-D248:133:H9FKKADXX:1
MP 80 HWI-D248:133:H9FKKADXX:1
MP 83 HISEQ:190:H9A78ADXX:1
MP 100 HWI-D248:133:H9FKKADXX:2
MP 104 HWI-D248:133:H9FKKADXX:2
MP 128 HISEQ:186:H9C50ADXX:1        
MP 218 HWI-D248:138:H9A7JADXX:2


CP 50 HWI-D248:132:H8LHJADXX:1
CP 51 HWI-D248:132:H8LHJADXX:2
CP 57 HWI-D248:132:H8LHJADXX:1
CP 58 HWI-D248:132:H8LHJADXX:1
CP 67 HISEQ:192:H9C0LADXX:1       
CP 68 HWI-D248:132:H8LHJADXX:1
CP 81 HWI-D248:133:H9FKKADXX:1
CP 82 HWI-D248:133:H9FKKADXX:1  
CP 84 HISEQ:190:H9A78ADXX:1
CP 98 HWI-D248:133:H9FKKADXX:2
CP 99 HWI-D248:133:H9FKKADXX:2



#To add
#RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample
#RGID=flowcell id:lane id RGLB=every one it’s own! RGPL=illumina RGPU={FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE} RGSM=SAMPLE
#Example #MP 45: RGID=H8LHJADXX:2 RGLB=MP45 RGPL=illumina RGPU=H8LHJADXX:2:45 RGSM=45


#45#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/45/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/45/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:2 RGLB=MP45 RGPL=illumina RGPU=H8LHJADXX:2:45 RGSM=45") | sbatch 
#47#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/47/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/47/rg_added_sorted.Aligned.out.bam RGID=H9C0LADXX:1 RGLB=MP47 RGPL=illumina RGPU=H9C0LADXX:1:47 RGSM=47") | sbatch  
#53#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/53/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/53/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:2 RGLB=MP53 RGPL=illumina RGPU=H8LHJADXX:2:53 RGSM=53") | sbatch 
#70#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/70/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/70/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:2 RGLB=MP70 RGPL=illumina RGPU=H8LHJADXX:2:70 RGSM=70") | sbatch
#71#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/71A/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/71A/rg_added_sorted.Aligned.out.bam RGID=H9C0LADXX:1 RGLB=MP71A RGPL=illumina RGPU=H9C0LADXX:1:71A RGSM=71A") | sbatch  
#78#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/78/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/78/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:1 RGLB=MP78 RGPL=illumina RGPU=H9FKKADXX:1:78 RGSM=78") | sbatch
#80#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/80/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/80/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:1 RGLB=MP80 RGPL=illumina RGPU=H9FKKADXX:1:80 RGSM=80") | sbatch
#83#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/83/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/83/rg_added_sorted.Aligned.out.bam RGID=H9A78ADXX:1 RGLB=MP83 RGPL=illumina RGPU=H9A78ADXX:1:83 RGSM=83") | sbatch
#100#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/100/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/100/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:2 RGLB=MP100 RGPL=illumina RGPU=H9FKKADXX:2:100 RGSM=100") | sbatch
#104#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/104/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/104/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:2 RGLB=MP104 RGPL=illumina RGPU=H9FKKADXX:2:104 RGSM=104") | sbatch
#128#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/128/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/128/rg_added_sorted.Aligned.out.bam RGID=H9C50ADXX:1 RGLB=MP128 RGPL=illumina RGPU=H9C50ADXX:1:128 RGSM=128") | sbatch
#218#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/218/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/218/rg_added_sorted.Aligned.out.bam RGID=H9A7JADXX:2 RGLB=MP218 RGPL=illumina RGPU=H9A7JADXX:2:218 RGSM=218") | sbatch


#50#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/50/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/50/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:1 RGLB=CP50 RGPL=illumina RGPU=H8LHJADXX:1:50 RGSM=50") | sbatch 
#51#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/51/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/51/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:2 RGLB=CP51 RGPL=illumina RGPU=H8LHJADXX:2:51 RGSM=51") | sbatch 
#57#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/57/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/57/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:1 RGLB=CP57 RGPL=illumina RGPU=H8LHJADXX:1:57 RGSM=57") | sbatch
#58#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/58/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/58/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:1 RGLB=CP58 RGPL=illumina RGPU=H8LHJADXX:1:58 RGSM=58") | sbatch
#67#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/67/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/67/rg_added_sorted.Aligned.out.bam RGID=H9C0LADXX:1  RGLB=CP67 RGPL=illumina RGPU=H9C0LADXX:1:67 RGSM=67") | sbatch
#68#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/68/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/68/rg_added_sorted.Aligned.out.bam RGID=H8LHJADXX:1 RGLB=CP68 RGPL=illumina RGPU=H8LHJADXX:1:68 RGSM=68") | sbatch
#81#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/81/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/81/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:1 RGLB=CP81 RGPL=illumina RGPU=H9FKKADXX:1:81 RGSM=81") | sbatch
#82#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/82/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/82/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:1 RGLB=CP82 RGPL=illumina RGPU=H9FKKADXX:1:82 RGSM=82") | sbatch
#84#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/84/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/84/rg_added_sorted.Aligned.out.bam RGID=H9A78ADXX:1 RGLB=CP84 RGPL=illumina RGPU=H9A78ADXX:1:84 RGSM=84") | sbatch
#98#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/98/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/98/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:2 RGLB=CP98 RGPL=illumina RGPU=H9FKKADXX:2:98 RGSM=98") | sbatch
#99#
(echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 1'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar AddOrReplaceReadGroups I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/99/sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/99/rg_added_sorted.Aligned.out.bam RGID=H9FKKADXX:2 RGLB=CP99 RGPL=illumina RGPU=H9FKKADXX:2:99 RGSM=99") | sbatch


#########MARK DUPLICATES (and make indexes)   
#Melpomene Adults    
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/
  mkdir ./$i/ 
  (echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 4'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar MarkDuplicates I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/$i/rg_added_sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/$i/dedupl.sorted.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=/data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/MP/$i/output.metrics") | sbatch
 done
#Cydno Adults
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/
  mkdir ./$i/ 
  (echo '#!/bin/bash'; echo '#SBATCH -J Picard'; echo '#SBATCH -n 4'; echo "java -jar /data/home/wolfproj/wolfproj-06/4_Picard_mappingQC/picard.jar MarkDuplicates I=/data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/$i/rg_added_sorted.Aligned.out.bam O=/data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/$i/dedupl.sorted.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=/data/home/wolfproj/wolfproj-06/6_Picard_dedupl/Adults/CP/$i/output.metrics") | sbatch
 done


#and so forth for the other samples
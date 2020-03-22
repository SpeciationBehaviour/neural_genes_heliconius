#Reannotate genomes using RABT (reference annotation based transcript), option -g #http://cole-trapnell-lab.github.io/cufflinks/cufflinks/


## Reannotate genome using transcripts from Adults

#Using melpomene RNA-reads
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/$i
  cd /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J cufflinks'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load cufflinks/2.1.1'; echo "cufflinks -g /data/home/wolfproj/wolfproj-06/Genome_annotations/Hmel2.5.gff3 /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/MPtoMP/$i/sorted.Aligned.out.bam -F 0.5 -I 20000 -A 0.4 –3-overhang-tolerance 5 –overlap-radius 5 –trim-3-dropoff-frac 0.3 –max-multiread-fraction 0.1 –min-frags-per-transfrag 50") | sbatch
 done

#COMMENTS
#-I #max intron length
#-A #Spliced reads with less than this percent of their length on each side of the junction are considered suspicious and are candidates for filtering prior to assembly. Default: 0.09.
#–max-bundle-length #default is 3,500,000 bp
#–trim-3-dropoff-frac The fraction of average coverage below which to trim the 3’ end of an assembled transcript. The default is 0.1.
#–3-overhang-tolerance The number of bp allowed to overhang the 3’ end of a reference transcript when determining if an assembled transcript should be merged with it (ie, the assembled transcript is not novel). The default is 600 bp.
#–overlap-radius Transfrags that are separated by less than this distance get merged together, and the gap is filled. Default: 50 (in bp).
#–min-frags-per-transfrag Assembled transfrags supported by fewer than this many aligned RNA-Seq fragments are not reported. (Default 10)
#-F #lower bound proportion isoform abundance

#Using cydno RNA-reads
cd /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/
individuals=$(ls -d *)    
for i in $individuals                 
 do
  mkdir /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/$i
  cd /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/$i
  (echo '#!/bin/bash'; echo '#SBATCH -J cufflinks'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load cufflinks/2.1.1'; echo "cufflinks -g /data/home/wolfproj/wolfproj-06/Genome_annotations/Hmel2.5.gff3 /data/home/wolfproj/wolfproj-06/3_STAR/Adults/2ndPass/CPtoMP/$i/sorted.Aligned.out.bam -F 0.5 -I 20000 -A 0.4 –3-overhang-tolerance 5 –overlap-radius 5 –trim-3-dropoff-frac 0.3 –max-multiread-fraction 0.1 –min-frags-per-transfrag 50") | sbatch
 done

# Combine annotations (gtf files) created individually into one single annotation

#first prepare file with list of files to combine
cd /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/combined_paramCufflink
nano combined_listparamCufflink.txt
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/100/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/104/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/128/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/218/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/45/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/47/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/53/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/70/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/71A/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/78/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/80/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_MP_paramCufflink/83/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/50/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/51/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/57/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/58/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/67/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/68/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/81/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/82/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/84/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/98/transcripts.gtf
/data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/ADULTS_CP_paramCufflink/99/transcripts.gtf

#Merge annotations
cd /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/combined_paramCufflink
(echo '#!/bin/bash'; echo '#SBATCH -J cuffmerge'; echo '#SBATCH -n 1'; echo '#SBATCH -t 23:59:59'; echo 'module load cufflinks/2.1.1'; echo "cuffmerge -g /data/home/wolfproj/wolfproj-06/Genome_annotations/Hmel2.5.gff3 -s /data/home/wolfproj/wolfproj-06/Genome_assemblies/Melpomene/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -o /data/home/wolfproj/wolfproj-06/14_guided_transcr_annotation/CUfflinks/combined_paramCufflink/combined_paramCufflink_stats combined_listparamCufflink.txt") | sbatch

#do the same for the other stages

#Run TrimGalore! (Wrapper script of Cutadapt) on RNA-seq reads, to trim low quality and adaptor bases + re-run FastQC (quality check) after having trimmed)

#example on  #Melpomene Adults
cd /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP
  mkdir ./$i/                 
  (echo '#!/bin/bash'; echo '#SBATCH -J Trim'; echo '#SBATCH -n 1'; echo 'module load cutadapt/1.9.1'; echo 'module load trimgalore/0.4.4'; echo 'module load fastqc/0.11.5'; echo "trim_galore --illumina --paired --length 30 --stringency 3 --trim-n --fastqc --output_dir /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/MP/$i/ /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/MP/$i/*_1.sanfastq.gz /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/MP/$i/*_2.sanfastq.gz") | sbatch
 done

####NOTES
#last two files R1 and R2 are the "paired-end" files
#note that illumina adapters (5’-3’) have the same first 13bp at both 5’ end (left of insert) and 3’end (right of insert), and TrimGalore use these first 13bp of the adapter to trim 
#quality cutoff is 20 (Phred score) by default (when trimming the 3’ end) #(cutadapt does only 3’ quality trimming) 
#--fastqc = Run FastQC once trimming is complete
#--length 30 (default =20 ) = both reads of a read-pair need to be longer than “length = bp” to be validated paired-end in the output file
#--stringency 3 # set at 3 (Cutadapt default)
#--trim-n = Removes Ns from either side of the read

#other example on #Cydno Adults
cd /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP
  mkdir ./$i/                 
  (echo '#!/bin/bash'; echo '#SBATCH -J Trim'; echo '#SBATCH -n 1'; echo 'module load cutadapt/1.9.1'; echo 'module load trimgalore/0.4.4'; echo 'module load fastqc/0.11.5'; echo "trim_galore --illumina --paired --length 30 --stringency 3 --trim-n --fastqc --output_dir /data/home/wolfproj/wolfproj-06/2_TrimGalore/results/Adults/CP/$i/ /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/CP/$i/*_1.sanfastq.gz /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/CP/$i/*_2.sanfastq.gz") | sbatch
 done

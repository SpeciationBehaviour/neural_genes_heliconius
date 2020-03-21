###send a fastqc (quality check) slurm job for every individual

#for Melpomene Adults
cd /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/MP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/1_FastQC/results/Adults/MP
  mkdir ./$i/                 
  (echo '#!/bin/bash'; echo '#SBATCH -J FastQC'; echo '#SBATCH -n 1'; echo 'module load fastqc/0.11.5'; echo "fastqc -o /data/home/wolfproj/wolfproj-06/1_FastQC/results/Adults/MP/$i/ -f fastq /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/Adults/MP/$i/*.sanfastq.gz") | sbatch
 done


## COMMENTS
## SBATCH -n 1 = 1 core
## last echo use "" otherwise treats $i as “character”
## LEGEND FastQC:  -o = output directory, -f = format, *.sanfastq.gz (at the end) are the files to check with fastqc 



#other example
#For Cydno 156h APF individuals
cd /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/156h_APF/CP
individuals=$(ls -d *)    
for i in $individuals                 
 do
  cd /data/home/wolfproj/wolfproj-06/1_FastQC/results/156h_APF/CP
  mkdir ./$i/                 
  (echo '#!/bin/bash'; echo '#SBATCH -J FastQC'; echo '#SBATCH -n 1'; echo 'module load fastqc/0.11.5'; echo "fastqc -o /data/home/wolfproj/wolfproj-06/1_FastQC/results/156h_APF/CP/$i/ -f fastq /data/home/wolfproj/wolfproj-06/GenePool_RNAseq_brain_data_July2014_Raw_reads/156h_APF/CP/$i/*.sanfastq.gz") | sbatch
 done


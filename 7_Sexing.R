#this script is for estimating HETEROZIGOSITY on the Z chromosome, in order to sex (pupal) samples

setwd("/Users/matteo/Desktop/")

#source("https://bioconductor.org/biocLite.R")
#biocLite("snpStats")
#biocLite("VariantAnnotation")
library(VariantAnnotation)
library(snpStats)

###### EXAMPLE WITH ADULT SAMPLES 

#CONVERT FILES with tabix/prodocue indexes in order to process them (on Cluster):
#module load samtools/1.4.1
#1) bgzip files
#cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/
#individuals=$(ls -d *)    
#for i in $individuals                 
 #do
  #cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/
  #bgzip -c $i.passed.output.vcf > $i.passed.output.vcf.gz 
 #done
#2) tabix
#cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/
#individuals=$(ls -d *)    
#for i in $individuals                 
 #do
  #cd /data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/$i/
  #tabix -p vcf $i.passed.output.vcf.gz  
#done

#DOWNLOAD FILES ON LAPTOP
#scp -r #LINK:/data/home/wolfproj/wolfproj-06/7_GATK/Adults/MP/*/*passed.output.vcf.gz.tbi /Users/matteo/Desktop/TBIs/Adults/MP/

#READ FILES
tab45 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/45.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/45.passed.output.vcf.gz.tbi")
tab47 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/47.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/47.passed.output.vcf.gz.tbi")
tab53 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/53.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/53.passed.output.vcf.gz.tbi")
tab70 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/70.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/70.passed.output.vcf.gz.tbi")
tab71 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/71A.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/71A.passed.output.vcf.gz.tbi")
tab78 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/78.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/78.passed.output.vcf.gz.tbi")
tab80 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/80.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/80.passed.output.vcf.gz.tbi")
tab83 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/83.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/83.passed.output.vcf.gz.tbi")
tab100 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/100.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/100.passed.output.vcf.gz.tbi")
tab104 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/104.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/104.passed.output.vcf.gz.tbi")
tab128 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/128.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/128.passed.output.vcf.gz.tbi")
tab218 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/MP/218.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/MP/218.passed.output.vcf.gz.tbi")

tab50 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/50.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/50.passed.output.vcf.gz.tbi")
tab51 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/51.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/51.passed.output.vcf.gz.tbi")
tab57 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/57.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/57.passed.output.vcf.gz.tbi")
tab58 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/58.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/58.passed.output.vcf.gz.tbi")
tab67 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/67.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/67.passed.output.vcf.gz.tbi")
tab68 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/68.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/68.passed.output.vcf.gz.tbi")
tab81 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/81.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/81.passed.output.vcf.gz.tbi")
tab82 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/82.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/82.passed.output.vcf.gz.tbi")
tab84 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/84.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/84.passed.output.vcf.gz.tbi")
tab98 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/98.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/98.passed.output.vcf.gz.tbi")
tab99 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/CP/99.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/CP/99.passed.output.vcf.gz.tbi")

#set the Z chromosome as the range for estimating heterozigosity
gr <- GRanges(seqnames = "Hmel221001o:1-13359691:+")   

vcf45 <- readVcf(tab45, "Melpomene",param=gr)
vcf47 <- readVcf(tab47, "Melpomene",param=gr)
vcf53 <- readVcf(tab53, "Melpomene",param=gr)
vcf70 <- readVcf(tab70, "Melpomene",param=gr)
vcf71 <- readVcf(tab71, "Melpomene",param=gr)
vcf78 <- readVcf(tab78, "Melpomene",param=gr)
vcf80 <- readVcf(tab80, "Melpomene",param=gr)
vcf83 <- readVcf(tab83, "Melpomene",param=gr)
vcf100 <- readVcf(tab100, "Melpomene",param=gr)
vcf104<- readVcf(tab104, "Melpomene",param=gr)
vcf128 <- readVcf(tab128, "Melpomene",param=gr)
vcf218 <- readVcf(tab218, "Melpomene",param=gr)

vcf50 <- readVcf(tab50, "Melpomene",param=gr)
vcf51 <- readVcf(tab51, "Melpomene",param=gr)
vcf57 <- readVcf(tab57, "Melpomene",param=gr)
vcf58 <- readVcf(tab58, "Melpomene",param=gr)
vcf67 <- readVcf(tab67, "Melpomene",param=gr)
vcf68 <- readVcf(tab68, "Melpomene",param=gr)
vcf81 <- readVcf(tab81, "Melpomene",param=gr)
vcf82 <- readVcf(tab82, "Melpomene",param=gr)
vcf84 <- readVcf(tab84, "Melpomene",param=gr)
vcf98 <- readVcf(tab98, "Melpomene",param=gr)
vcf99 <- readVcf(tab99, "Melpomene",param=gr)

#exlcude all but biallelic SNPs 
snpmat45 <- genotypeToSnpMatrix(vcf45)    
snpmat47 <- genotypeToSnpMatrix(vcf47)
snpmat53 <- genotypeToSnpMatrix(vcf53)
snpmat70 <- genotypeToSnpMatrix(vcf70)
snpmat71 <- genotypeToSnpMatrix(vcf71)
snpmat78 <- genotypeToSnpMatrix(vcf78)
snpmat80 <- genotypeToSnpMatrix(vcf80)
snpmat83 <- genotypeToSnpMatrix(vcf83)
snpmat100 <- genotypeToSnpMatrix(vcf100)
snpmat104 <- genotypeToSnpMatrix(vcf104)
snpmat128 <- genotypeToSnpMatrix(vcf128)
snpmat218 <- genotypeToSnpMatrix(vcf218)

snpmat50 <- genotypeToSnpMatrix(vcf50)
snpmat51 <- genotypeToSnpMatrix(vcf51)
snpmat57 <- genotypeToSnpMatrix(vcf57)
snpmat58 <- genotypeToSnpMatrix(vcf58)
snpmat67 <- genotypeToSnpMatrix(vcf67)
snpmat68 <- genotypeToSnpMatrix(vcf68)
snpmat81 <- genotypeToSnpMatrix(vcf81)
snpmat82 <- genotypeToSnpMatrix(vcf82)
snpmat84 <- genotypeToSnpMatrix(vcf84)
snpmat98 <- genotypeToSnpMatrix(vcf98)
snpmat99 <- genotypeToSnpMatrix(vcf99)

#estimate heterozigosity
summary(snpmat45$genotypes)  #Heterozygosity #0.4872 #0.49     #45 M
summary(snpmat47$genotypes)   #0.457 #0.46                   #47 M
summary(snpmat53$genotypes)    #0.04281  #0.04               #53 F
summary(snpmat70$genotypes)    #0.4684 #0.47              #70 M
summary(snpmat71$genotypes)   #0.4782 # 0.48              #71 M          
summary(snpmat78$genotypes)   #0.04745 # 0.05             #78 F
summary(snpmat80$genotypes)  #0.03883 #0.04             #80 F
summary(snpmat83$genotypes)   #0.4631 #0.46            #83 M
summary(snpmat100$genotypes)   #0.486 #0.49              #100 M
summary(snpmat104$genotypes)    #0.4624 #0.46           #104 M
summary(snpmat128$genotypes)    #0.04477 #0.04           #128 F
summary(snpmat218$genotypes)    #0.04637 # 0.05             #218 F

summary(snpmat50$genotypes)  #0.02303 #0.02      #50 F
summary(snpmat51$genotypes)   #0.02328 #0.02      #51 F
summary(snpmat57$genotypes)    #0.2331 #0.23      #57 M
summary(snpmat58$genotypes)    #0.02303 #0.02     #58 F   
summary(snpmat67$genotypes)   #0.01585 #0.02      #67 F         
summary(snpmat68$genotypes)   #0.01527 #0.01      #68 F
summary(snpmat81$genotypes)   #0.02415 #0.02      #81 F
summary(snpmat82$genotypes)   #0.2325  #0.23    #82 M
summary(snpmat84$genotypes)    #0.02071 #0.02      #84 F
summary(snpmat98$genotypes)     #0.2539 #0.25    #98 M
summary(snpmat99$genotypes)    #0.2438  #0.24    #99 M



####### Other example with F1 hybrids files
#read files
tab42 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/F1/42.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/F1/42.passed.output.vcf.gz.tbi")
tab49 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/F1/49.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/F1/49.passed.output.vcf.gz.tbi")
tab56 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/F1/56.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/F1/56.passed.output.vcf.gz.tbi")
tab69 <- TabixFile("/Users/matteo/Desktop/GZs/Adults/F1/69.passed.output.vcf.gz","/Users/matteo/Desktop/TBIs/Adults/F1/69.passed.output.vcf.gz.tbi")
gr <- GRanges(seqnames = "Hmel221001o:1-13359691:+")   #set Z chr as the range
vcf42 <- readVcf(tab42, "Melpomene",param=gr)
vcf49 <- readVcf(tab49, "Melpomene",param=gr)
vcf56 <- readVcf(tab56, "Melpomene",param=gr)
vcf69 <- readVcf(tab69, "Melpomene",param=gr)
snpmat42 <- genotypeToSnpMatrix(vcf42)    #restrict to biallelic SNPs
snpmat49 <- genotypeToSnpMatrix(vcf49)
snpmat56 <- genotypeToSnpMatrix(vcf56)
snpmat69 <- genotypeToSnpMatrix(vcf69)
#estimate heterozigosity
summary(snpmat42$genotypes)  #0.6378 #0.64    #42 M
summary(snpmat49$genotypes)  #0.6358 #0.64    #49 M
summary(snpmat56$genotypes)  #0.05236 #0.05     #56 F
summary(snpmat69$genotypes)  #0.04294 #0.04    #69 F

#and so for all other samples

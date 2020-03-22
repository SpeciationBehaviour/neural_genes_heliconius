#retrieve predicted protein sequences of the Hmel2.5 genome at http://butterflygenome.org/?q=node/4
#modify file “Heliconius_melpomene_melpomene_Hmel2.5.proteins.fa” to remove all * (stop codons)
sed '/^>/ s/ .*//' Heliconius_melpomene_melpomene_Hmel2.5.proteins.fa > Heliconius_melpomene_melpomene_Hmel2.5.proteins.fasta

#download InterProScan
##run Interproscan on protein sequences
cd /data/home/wolfproj/wolfproj-06/8_InterProScan/
(echo '#!/bin/bash'; echo '#SBATCH -J interproscan'; echo '#SBATCH -n 5'; echo "/data/home/wolfproj/wolfproj-06/8_InterProScan/interproscan-5.27-66.0/interproscan.sh -i /data/home/wolfproj/wolfproj-06/8_InterProScan/Heliconius_melpomene_melpomene_Hmel2.5.proteins.fasta -iprlookup -goterms --pathways -b /data/home/wolfproj/wolfproj-06/8_InterProScan/RESULTS/") | sbatch  
#Databases used:
#Coils-2.2.1,Gene3D-4.1.0,Hamap-2017_10,MobiDBLite-1.0,Pfam-31.0,PIRSF-3.02,PRINTS-42.0,ProDom-2006.1,ProSitePatterns-2017_09,ProSiteProfiles-2017_09,SFLD-3,SMART-7.1,SUPERFAMILY-1.75,TIGRFAM-15.0

#Download PANTHER
##run Panther on protein sequences
cd /data/home/wolfproj/wolfproj-06/8_InterProScan/
(echo '#!/bin/bash'; echo '#SBATCH -J panther'; echo '#SBATCH -n 2'; echo "/data/home/wolfproj/wolfproj-06/8_InterProScan/interproscan-5.27-66.0/interproscan.sh -i /data/home/wolfproj/wolfproj-06/8_InterProScan/Heliconius_melpomene_melpomene_Hmel2.5.proteins.fasta -appl PANTHER -iprlookup -goterms --pathways -b /data/home/wolfproj/wolfproj-06/8_InterProScan/PANTHER_RESULTS/") | sbatch

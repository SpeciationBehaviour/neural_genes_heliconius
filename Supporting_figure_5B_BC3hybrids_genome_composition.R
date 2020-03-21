setwd("/Users/matteo/Desktop/")




######################################    SNPs density

#Cydno "Pure", 156h
#4
MPcoverage_4 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/4.MP_coverage.txt",header=FALSE)
CPcoverage_4 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/4.CP_coverage.txt",header=FALSE)
#rename column
snps4<-cbind(MPcoverage_4,CPcoverage_4)
colnames(snps4)[8] <- "V5"

#F1 hybrids
#42 (male)
MPcoverage_42 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/42.MP_coverage.txt",header=FALSE)
CPcoverage_42 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/42.CP_coverage.txt",header=FALSE)
#
snps42<-cbind(MPcoverage_42,CPcoverage_42)
colnames(snps42)[8] <- "V5"


#156h hybrids
MPcoverage_116 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/116.MP_coverage.txt",header=FALSE)
CPcoverage_116 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/116.CP_coverage.txt",header=FALSE)
MPcoverage_123 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/123.MP_coverage.txt",header=FALSE)
CPcoverage_123 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/123.CP_coverage.txt",header=FALSE)
MPcoverage_126 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/126.MP_coverage.txt",header=FALSE)
CPcoverage_126 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/126.CP_coverage.txt",header=FALSE)
MPcoverage_131 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/131.MP_coverage.txt",header=FALSE)
CPcoverage_131 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/131.CP_coverage.txt",header=FALSE)
MPcoverage_139 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/139.MP_coverage.txt",header=FALSE)
CPcoverage_139 <-read.table("/Users/matteo/Desktop/snps_density_filtered_RNA/139.CP_coverage.txt",header=FALSE)

#
snps116<-cbind(MPcoverage_116,CPcoverage_116)
colnames(snps116)[8] <- "V5"
#
snps123<-cbind(MPcoverage_123,CPcoverage_123)
colnames(snps123)[8] <- "V5"
#
snps126<-cbind(MPcoverage_126,CPcoverage_126)
colnames(snps126)[8] <- "V5"
#
snps131<-cbind(MPcoverage_131,CPcoverage_131)
colnames(snps131)[8] <- "V5"
#
snps139<-cbind(MPcoverage_139,CPcoverage_139)
colnames(snps139)[8] <- "V5"



######################################################################################################PLOTS
###########################################################################
library(ggplot2)
library(gridExtra)

#AREA PLOT #with FRACTIONS
#
colnames(snps4)[5] <- "V5"
colnames(snps4)[6] <- "V6"
colnames(snps4)[7] <- "V7"
colnames(snps4)[8] <- "V8"
#
colnames(snps42)[5] <- "V5"
colnames(snps42)[6] <- "V6"
colnames(snps42)[7] <- "V7"
colnames(snps42)[8] <- "V8"
#
colnames(snps116)[5] <- "V5"
colnames(snps116)[6] <- "V6"
colnames(snps116)[7] <- "V7"
colnames(snps116)[8] <- "V8"
#
colnames(snps123)[5] <- "V5"
colnames(snps123)[6] <- "V6"
colnames(snps123)[7] <- "V7"
colnames(snps123)[8] <- "V8"
#
colnames(snps126)[5] <- "V5"
colnames(snps126)[6] <- "V6"
colnames(snps126)[7] <- "V7"
colnames(snps126)[8] <- "V8"
#
colnames(snps131)[5] <- "V5"
colnames(snps131)[6] <- "V6"
colnames(snps131)[7] <- "V7"
colnames(snps131)[8] <- "V8"
#
colnames(snps139)[5] <- "V5"
colnames(snps139)[6] <- "V6"
colnames(snps139)[7] <- "V7"
colnames(snps139)[8] <- "V8"



###########add fractions column
#
snps4$V9 <- snps4$V4 / (snps4$V4 + snps4$V8)
snps4$V10 <- snps4$V8 / (snps4$V4 + snps4$V8)
#
snps42$V9 <- snps42$V4 / (snps42$V4 + snps42$V8)
snps42$V10 <- snps42$V8 / (snps42$V4 + snps42$V8)

#
snps116$V9 <- snps116$V4 / (snps116$V4 + snps116$V8)
snps116$V10 <- snps116$V8 / (snps116$V4 + snps116$V8)
#
snps123$V9 <- snps123$V4 / (snps123$V4 + snps123$V8)
snps123$V10 <- snps123$V8 / (snps123$V4 + snps123$V8)
#
snps126$V9 <- snps126$V4 / (snps126$V4 + snps126$V8)
snps126$V10 <- snps126$V8 / (snps126$V4 + snps126$V8)
#
snps131$V9 <- snps131$V4 / (snps131$V4 + snps131$V8)
snps131$V10 <- snps131$V8 / (snps131$V4 + snps131$V8)
#
snps139$V9 <- snps139$V4 / (snps139$V4 + snps139$V8)
snps139$V10 <- snps139$V8 / (snps139$V4 + snps139$V8)


#######add column sum of snps
snps4$control<-snps4$V4 + snps4$V8
snps42$control<-snps42$V4 + snps42$V8

snps116$control<-snps116$V4 + snps116$V8
snps123$control<-snps123$V4 + snps123$V8
snps126$control<-snps126$V4 + snps126$V8
snps131$control<-snps131$V4 + snps131$V8
snps139$control<-snps139$V4 + snps139$V8




#remove values  where less than 30 snps 
snps4<-snps4[snps4$control>30,] 
snps42<-snps42[snps42$control>30,] 

snps116<-snps116[snps116$control>30,] 
snps123<-snps123[snps123$control>30,] 
snps126<-snps126[snps126$control>30,] 
snps131<-snps131[snps131$control>30,] 
snps139<-snps139[snps139$control>30,] 



#order chromosomes
snps4$chr_f = factor(snps4$V1, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21'))
snps42$chr_f = factor(snps42$V1, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21'))

snps116$chr_f = factor(snps116$V1, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21'))
snps123$chr_f = factor(snps123$V1, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21'))
snps126$chr_f = factor(snps126$V1, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21'))
snps131$chr_f = factor(snps131$V1, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21'))
snps139$chr_f = factor(snps139$V1, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21'))




#PLOTS
#
ggplot(snps4, aes(x=V2, y=V9)) + geom_area(color="darkgoldenrod1",fill = "darkgoldenrod1",alpha = 0.5) +
  geom_area(data=snps4, aes(x=V2, y=V10),color="cadetblue3",fill = "cadetblue3",alpha = 0.5) +  theme_void() + facet_grid(~chr_f) + theme(title=element_blank()) 
#
ggplot(snps42, aes(x=V2, y=V9)) + geom_area(color="darkgoldenrod1",fill = "darkgoldenrod1",alpha = 0.5) +
  geom_area(data=snps42, aes(x=V2, y=V10),color="cadetblue3",fill = "cadetblue3",alpha = 0.5) +  theme_void() + facet_grid(~chr_f) + theme(title=element_blank()) 

#156h hybrids subset
#
ggplot(snps116, aes(x=V2, y=V9)) + geom_area(color="darkgoldenrod1",fill = "darkgoldenrod1",alpha = 0.5) +
  geom_area(data=snps116, aes(x=V2, y=V10),color="cadetblue3",fill = "cadetblue3",alpha = 0.5) +  theme_void() + facet_grid(~chr_f) + theme(title=element_blank()) 
#
ggplot(snps123, aes(x=V2, y=V9)) + geom_area(color="darkgoldenrod1",fill = "darkgoldenrod1",alpha = 0.5) +
  geom_area(data=snps123, aes(x=V2, y=V10),color="cadetblue3",fill = "cadetblue3",alpha = 0.5) +  theme_void() + facet_grid(~chr_f) + theme(title=element_blank())
#
ggplot(snps126, aes(x=V2, y=V9)) + geom_area(color="darkgoldenrod1",fill = "darkgoldenrod1",alpha = 0.5) +
  geom_area(data=snps126, aes(x=V2, y=V10),color="cadetblue3",fill = "cadetblue3",alpha = 0.5) +  theme_void() + facet_grid(~chr_f) + theme(title=element_blank()) 
#
ggplot(snps131, aes(x=V2, y=V9)) + geom_area(color="darkgoldenrod1",fill = "darkgoldenrod1",alpha = 0.5) +
  geom_area(data=snps131, aes(x=V2, y=V10),color="cadetblue3",fill = "cadetblue3",alpha = 0.5) +  theme_void() + facet_grid(~chr_f) + theme(title=element_blank())
#
ggplot(snps139, aes(x=V2, y=V9)) + geom_area(color="darkgoldenrod1",fill = "darkgoldenrod1",alpha = 0.5) +
  geom_area(data=snps139, aes(x=V2, y=V10),color="cadetblue3",fill = "cadetblue3",alpha = 0.5) +  theme_void() + facet_grid(~chr_f) + theme(title=element_blank())









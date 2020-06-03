setwd("/Users/matteo/Desktop/tables_allele_specific/toCP")
#setwd("/Users/matteo/Desktop/tables_allele_specific/toMP") #no informative spns for candidates

library(dplyr)
library(GenomicRanges)
library(DESeq2)


####read tables produced with previous script (with number of reads in F1 hybrids mapping to the cydno (ALT) or melpomene allele (REF)
F1_42<-read.csv("42.allele_counts.csv",sep="\t",header=TRUE) 
F1_49<-read.csv("49.allele_counts.csv",sep="\t",header=TRUE) 
F1_56<-read.csv("56.allele_counts.csv",sep="\t",header=TRUE) 
F1_69<-read.csv("69.allele_counts.csv",sep="\t",header=TRUE) 

F1_90<-read.csv("90.allele_counts.csv",sep="\t",header=TRUE)
F1_152<-read.csv("152.allele_counts.csv",sep="\t",header=TRUE)
F1_200<-read.csv("200.allele_counts.csv",sep="\t",header=TRUE)
F1_273<-read.csv("273.allele_counts.csv",sep="\t",header=TRUE)
F1_292<-read.csv("292.allele_counts.csv",sep="\t",header=TRUE)
F1_298<-read.csv("298.allele_counts.csv",sep="\t",header=TRUE)


#keep onbly certain columns 1,2,6,7,8
F1_42<-F1_42[,c(1,2,6,7,8)]
F1_49<-F1_49[,c(1,2,6,7,8)]
F1_56<-F1_56[,c(1,2,6,7,8)]
F1_69<-F1_69[,c(1,2,6,7,8)]
#rename
colnames(F1_42)[1] <- "scaffold"
colnames(F1_49)[1] <- "scaffold"
colnames(F1_56)[1] <- "scaffold"
colnames(F1_69)[1] <- "scaffold"

colnames(F1_90)[1] <- "scaffold"
colnames(F1_152)[1] <- "scaffold"
colnames(F1_200)[1] <- "scaffold"
colnames(F1_273)[1] <- "scaffold"
colnames(F1_292)[1] <- "scaffold"
colnames(F1_298)[1] <- "scaffold"


#COMBINE FILES
new<-full_join(F1_42,F1_49, by=c("position","scaffold"))  
new2<-full_join(new,F1_56, by=c("position","scaffold"))
new3<-full_join(new2,F1_69, by=c("position","scaffold"))
new4<-full_join(new3,F1_90, by=c("position","scaffold"))
new5<-full_join(new4,F1_152, by=c("position","scaffold"))
new6<-full_join(new5,F1_200, by=c("position","scaffold"))
new7<-full_join(new6,F1_273, by=c("position","scaffold"))
new8<-full_join(new7,F1_292, by=c("position","scaffold"))
new9<-full_join(new8,F1_298, by=c("position","scaffold"))



#ASSIGN SNPs to ANNOTATED GENES
setwd("/Users/matteo/Desktop/")

annotation<-read.table("adults_Cufflinks_gene_info.txt",header=TRUE)     
#need to add chromosome information
GRgenes<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))

#change column names and add "end" position (same as start), needed by GRanges
colnames(new9)[1] <- "scaffold"
colnames(new9)[2] <- "start"
new9$end = new9$start
new<-new9[,c(1,2,81,3:4,6:7,9:10,12:13,18:19,29:30,40:41,51:52,62:63,73:74)]
GRsnps<-makeGRangesFromDataFrame(new, seqnames.field=c("scaffold"))

#find overalps between snps and genes
olaps <- findOverlaps(GRgenes, GRsnps)
a<-cbind(annotation[queryHits(olaps),], new[subjectHits(olaps),])


### SUM snps counts into allele counts 
#replace NAs
a[is.na(a)] <- 0
F1_42_ref<-rowsum(a$refCount.x, a$gene_id)
F1_42_alt<-rowsum(a$altCount.x, a$gene_id)
F1_49_ref<-rowsum(a$refCount.y, a$gene_id)
F1_49_alt<-rowsum(a$altCount.y, a$gene_id)
F1_56_ref<-rowsum(a$refCount.x.x, a$gene_id)
F1_56_alt<-rowsum(a$altCount.x.x, a$gene_id)
F1_69_ref<-rowsum(a$refCount.y.y, a$gene_id)
F1_69_alt<-rowsum(a$altCount.y.y, a$gene_id)
F1_90_ref<-rowsum(a$refCount.x.x.x, a$gene_id)
F1_90_alt<-rowsum(a$altCount.x.x.x, a$gene_id)
F1_152_ref<-rowsum(a$refCount.y.y.y, a$gene_id)
F1_152_alt<-rowsum(a$altCount.y.y.y, a$gene_id)
F1_200_ref<-rowsum(a$refCount.x.x.x.x, a$gene_id)
F1_200_alt<-rowsum(a$altCount.x.x.x.x, a$gene_id)
F1_273_ref<-rowsum(a$refCount.y.y.y.y, a$gene_id)
F1_273_alt<-rowsum(a$altCount.y.y.y.y, a$gene_id)
F1_292_ref<-rowsum(a$refCount.x.x.x.x.x, a$gene_id)
F1_292_alt<-rowsum(a$altCount.x.x.x.x.x, a$gene_id)
F1_298_ref<-rowsum(a$refCount.y.y.y.y.y, a$gene_id)
F1_298_alt<-rowsum(a$altCount.y.y.y.y.y, a$gene_id)


matrix<-cbind(F1_42_ref,F1_42_alt,F1_49_ref,F1_49_alt,F1_56_ref,F1_56_alt,F1_69_ref,F1_69_alt,
              F1_90_ref,F1_90_alt,F1_152_ref,F1_152_alt,F1_200_ref,F1_200_alt,
              F1_273_ref,F1_273_alt,F1_292_ref,F1_292_alt,F1_298_ref,F1_298_alt)

table<-data.frame(matrix, rownames = TRUE)

table <- add_rownames(table, "gene_id")
colnames(table)[2] <- "42_refCount"
colnames(table)[3] <- "42_altCount"
colnames(table)[4] <- "49_refCount"
colnames(table)[5] <- "49_altCount"
colnames(table)[6] <- "56_refCount"
colnames(table)[7] <- "56_altCount"
colnames(table)[8] <- "69_refCount"
colnames(table)[9] <- "69_altCount"
colnames(table)[10] <- "90_refCount"
colnames(table)[11] <- "90_altCount"
colnames(table)[12] <- "152_refCount"
colnames(table)[13] <- "152_altCount"
colnames(table)[14] <- "200_refCount"
colnames(table)[15] <- "200_altCount"
colnames(table)[16] <- "273_refCount"
colnames(table)[17] <- "273_altCount"
colnames(table)[18] <- "292_refCount"
colnames(table)[19] <- "292_altCount"
colnames(table)[20] <- "298_refCount"
colnames(table)[21] <- "298_altCount"


##note on Cufflinks transcript-based annotation
#CUFF.13201 = Ionotropic glut receptor
#CUFF.13212 = regucalcin2


##### CONDUCT TEST for DIFFERENTIAL ALLELE expression with DEseq2 approach (model= ~0+individual+allele)
table2<-table[,c(2,4,6,8,10,12,14,16,18,20,3,5,7,9,11,13,15,17,19,21)]
rownames(table2) <- table$gene_id
table2<- as.matrix(table2)
(individual <- factor(c(1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10)))
(allele <- factor(c(rep("mel", 10), rep("cydno", 10))))
(coldata <- data.frame(row.names=colnames(table2), individual,allele))
dds <- DESeqDataSetFromMatrix(countData=table2, colData=coldata, design=~0+individual+allele)
sizeFactors(dds) <- 1    #should not normalize between samples
dds <- DESeq(dds)
res <- results(dds)
#write.table(res, file="Alleles_newdata_ALL_expr.tx",  quote=FALSE)  #add gene_id file first column


### TO COMPARE TO SPECIES LEVEL
allAdult<-read.table("Adults_toCufflinks_genecounts.txt",header=TRUE)
allAdult<-allAdult[,c(1:24)]
MPnew<-read.table("merged_htseq_MPtoCuff_genecounts.txt",header=TRUE)
CPnew<-read.table("merged_htseq_CPtoCuff_genecounts.txt",header=TRUE)
all <- merge(allAdult,MPnew,by="gene_id")
all2 <- merge(all,CPnew,by="gene_id")
#write.table(all2, file="All_species_samples_counts.txt",  quote=FALSE, row.names = FALSE)  #add gene_id file first column

allAdult<-read.table("All_species_samples_counts.txt",header=TRUE, row.names = 1)
allAdult<- as.matrix(allAdult)
(sex <- factor(c(rep("male",2), rep("female", 1),rep("male",2), rep("female",2),rep("male",3),rep("female",4),rep("male",1),rep("female",4),rep("male",1),rep("female",1),rep("male",12))))
(batch <- factor(c(rep("2014",23), rep("2019", 10))))
(conditionAdult <- factor(c(rep("mel", 12), rep("cydno", 11),rep("mel",5),rep("cydno",5))))
(coldata <- data.frame(row.names=colnames(allAdult), conditionAdult,sex,batch))
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+batch+conditionAdult)
ddsSex <- DESeq(ddsSex)
res <- results(ddsSex)
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.05)
#write.table(res, file="species_data_ALL_expr.txt",  quote=FALSE)  #add gene_id file first column



###### TO MAKE FIGURE 3
#CUFF.13212  #regucalcin2
#hybrids levels
#value log2 fold-change #standard error of the (logarithmic) fold change (lfcSE)
#1.76218 #0.4288
#species level #(with random factor sequencing batch)
#3.4645 #0.290
#CUFF.13201 #Ionotropic glut receptor
#hybrids levels
#1.10593 #0.1907
#species level
#2.099 #0.1882
d = data.frame(
  x  = c(1.76218)
  , y  = c(3.4645)
  , sdx = c(0.4288)
  , sdy = c(0.290),
  x2  = c(1.10593)
  , y2  = c(2.099 )
  , sdx2 = c(0.1907)
  , sdy2 = c(0.1882)
)

#
plot (d$x, d$y, ylim=c(-4,4),xlim=c(-4,4),ylab="",xlab="")
# xo, yo, x1,y1
segments(d$x - d$sdx,d$y,d$x + d$sdx,d$y )
segments(d$x,d$y - d$sdy,d$x,d$y + d$sdy)
points(d$x2, d$y2)
segments(d$x2 - d$sdx2,d$y2,d$x2 + d$sdx2,d$y2 )
segments(d$x2,d$y2 - d$sdy2,d$x2,d$y2 + d$sdy2)

abline(0,0)
abline (v=0)
abline (v=1, lty=2,col="gray62")
abline (1,0, lty=2,col="gray62")
abline(0,1, lty=2,col="gray87")


### TO MAKE SUPPORTING FIGURE 7
hybrids<-read.table("Alleles__newdata_ALL_expr.txt",header=TRUE) 
hist(hybrids$log2FoldChange, xlim=c(-4,4),breaks=40, col = c(c(rep("gold", 24),c(rep("dodgerblue", 17)))))

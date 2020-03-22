setwd("/Users/matteo/Desktop/")

library(dplyr)
library(GenomicRanges)
library(DESeq2)

####read tables produced with previous script (with number of reads in F1 hybrids mapping to the cydno (ALT) or melpomene allele (REF)
F1_42<-read.csv("42.allele_counts.csv",sep="\t",header=TRUE) 
F1_49<-read.csv("49.allele_counts.csv",sep="\t",header=TRUE) 
F1_56<-read.csv("56.allele_counts.csv",sep="\t",header=TRUE) 
F1_69<-read.csv("69.allele_counts.csv",sep="\t",header=TRUE) 
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
#COMBINE FILES
new<-full_join(F1_42,F1_49, by=c("position","scaffold"))  
new2<-full_join(new,F1_56, by=c("position","scaffold"))
new3<-full_join(new2,F1_69, by=c("position","scaffold"))


#ASSIGN SNPs to ANNOTATED GENES
annotation<-read.table("adults_gene_info.txt",header=TRUE)     
#need to add chromosome information
GRgenes<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))

#change column names and add "end" position (same as start), needed by GRanges
colnames(new3)[1] <- "scaffold"
colnames(new3)[2] <- "start"
new3$end = new3$start
new3<-new3[,c(1,2,15,3:14)]
GRsnps<-makeGRangesFromDataFrame(new3, seqnames.field=c("scaffold"))

#find overalps between snps and genes
olaps <- findOverlaps(GRgenes, GRsnps)
a<-cbind(annotation[queryHits(olaps),], new3[subjectHits(olaps),])


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
matrix<-cbind(F1_42_ref,F1_42_alt,F1_49_ref,F1_49_alt,F1_56_ref,F1_56_alt,F1_69_ref,F1_69_alt)
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

##note on Cufflinks transcript-based annotation
#CUFF.13201 = Ionotropic glut receptor
#CUFF.13212 = regucalcin2


##### CONDUCT TEST for DIFFERENTIAL ALLELE expression with DEseq2 approach (model= ~0+individual+allele)
table2<-table[,c(1,2,4,6,8,3,5,7,9)]
rownames(table2) <- table$gene_id
table2$gene_id<-NULL
table2<- as.matrix(table2)
(individual <- factor(c(1,2,3,4,1,2,3,4)))
(allele <- factor(c(rep("mel", 4), rep("cydno", 4))))
(coldata <- data.frame(row.names=colnames(table2), individual,allele))
dds <- DESeqDataSetFromMatrix(countData=table2, colData=coldata, design=~0+individual+allele)
sizeFactors(dds) <- 1    #should not normalize between samples
dds <- DESeq(dds)
res <- results(dds)
#write.table(res, file="Alleles_withcluster_newannALL_expr.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column




######### TO MAKE FIGURE 3

#CUFF.13212  #regucalcin2
#hybrids levels
#value log2 fold-change #standard error of the (logarithmic) fold change (lfcSE)
#1.127865 # 1.383361
#species level
#3.826666 #0.3493866

#CUFF.13201 #Ionotropic glut receptor
#hybrids levels
#1.298138 #0.5475794
#species level
#2.21346 #0.2067517

#######
d = data.frame(
  x  = c(1.127865)
  , y  = c(3.826666)
  , sdx = c(1.383361)
  , sdy = c(0.3493866),
  x2  = c(1.298138)
  , y2  = c(2.21346)
  , sdx2 = c(0.5475794)
  , sdy2 = c(0.2067517)
)

#Ionotropic glut receptor
d = data.frame(
  x2  = c(1.298138)
  , y2  = c(2.21346)
  , sdx2 = c(0.5475794)
  , sdy2 = c(0.2067517)
)
#regucalcin2
d = data.frame(
  x  = c(1.127865)
  , y  = c(3.826666)
  , sdx = c(1.383361)
  , sdy = c(0.3493866))

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
#abline (v=1, lty=2,col="gold")
#abline (v=-1, lty=2,col="blue")
#abline (1,0, lty=2,col="gold")
#abline (-1,0, lty=2,col="blue")
abline(0,1)


######### TO MAKE SUPPORTING FIGURE 7
hybrids<-read.table("Alleles_withcluster_newannALL_expr.txt",header=TRUE) 
hist(hybrids$log2FoldChange, xlim=c(-4,4),breaks=40, col = c(c(rep("yellow", 20),c(rep("blue", 17)))))

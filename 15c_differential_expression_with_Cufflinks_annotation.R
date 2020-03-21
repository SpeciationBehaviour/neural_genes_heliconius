setwd("/Users/matteo/Desktop/")  #

library(DESeq2)
library(dplyr)


##############################       Adults

###new annotation
#extract from ANNOTATION gene and scaffold info #MELPOMENE
#read annotation #and make tables with gene, scaffold, start, end, INFO
#gff3<-read.delim("Adults_combined_Cufflinks_transcripts.gtf", header=F, comment.char="#")
#genesattribute <- data.frame(do.call('rbind', strsplit(as.character(gff3$V9),';',fixed=TRUE)))
#genesattribute2 <- data.frame(do.call('rbind', strsplit(as.character(genesattribute$X1),' ',fixed=TRUE)))
#gff3<-cbind(gff3,genesattribute2)
#annotation<-gff3[c(11,1,4,5)]
#colnames(annotation)[1] <- "gene_id"
#colnames(annotation)[2] <- "scaffold"
#colnames(annotation)[3] <- "start"
#colnames(annotation)[4] <- "end"
#genes<-annotation[!duplicated(annotation$gene_id),]
#write.table(genes, file="adults_Cufflinks_gene_info.txt", row.names = FALSE, quote=FALSE)  
annotation<-read.table("adults_Cufflinks_gene_info.txt",header=TRUE) 


#read output files from HTseq-count
#mapped to MP
MPtoMP <-read.table("merged_htseq_MPtoCuff_genecounts.txt",header=TRUE)
CPtoMP <-read.table("merged_htseq_CPtoCuff_genecounts.txt",header=TRUE)
F1toMP <-read.table("merged_htseq_F1toCuff_genecounts.txt",header=TRUE)
all <- merge(MPtoMP,CPtoMP,by="gene_id")
all <- merge(all,F1toMP,by="gene_id")
#write.table(all, file="Adults_toCufflinks_genecounts.txt", row.names = FALSE, quote=FALSE)



#Multi-factor design #sex as random factor   
allAdult<-read.table("Adults_toCufflinks_genecounts.txt",header=TRUE, row.names = 1)
allAdult<-allAdult[,c(1:23)]  
allAdult<- as.matrix(allAdult)
(sex <- factor(c(rep("male",2), rep("female", 1),rep("male",2), rep("female",2),rep("male",3), rep("female",4),rep("male",1), rep("female",4),rep("male",1), rep("female",1),rep("male",2))))
(conditionAdult <- factor(c(rep("mel", 12), rep("cydno", 11))))
(coldata <- data.frame(row.names=colnames(allAdult), conditionAdult,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+conditionAdult)
ddsSex <- DESeq(ddsSex)
res <- results(ddsSex)
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.05)
#write.table(resSig, file="Adults_DE_sexfactor_toCufflinks.txt", row.names = TRUE, quote=FALSE) 
#add gene_id file first column

#FOR ALL GENES 
#write.table(res, file="Adults_sexfactorALL_toCufflinks.txt", row.names = TRUE, quote=FALSE) 



diffAdult<-read.table("Adults_DE_sexfactor_toCufflinks.txt",header=TRUE)
upmel<-diffAdult[diffAdult$log2FoldChange >1,]  #
upcydno<-diffAdult[diffAdult$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")
#CHR18
only18_1<-all[all$scaffold == "Hmel218001o",] 
only18_1<-only18_1[only18_1$start >8997,] #...
only18_2<-all[all$scaffold == "Hmel218002o",] 
only18_3<-all[all$scaffold =="Hmel218003o",] #  
onlypeak<-only18_3[only18_3$end <275070,] #  
onlylodscore<-only18_3[only18_3$end <2393341,] 
#CUFF.13201 #IgluR
##CUFF.13212 #regucalcin


###########     F1 vs cydno
#reread allAdult
allAdult<-read.table("Adults_toCufflinks_genecounts.txt",header=TRUE, row.names = 1)
allAdult<-allAdult[,c(13:24,26:28)] 
allAdult<- as.matrix(allAdult)
(sex <- factor(c(rep("male",2), rep("female", 1),rep("male",2), rep("female",2),rep("male",3), rep("female",4),rep("male",1), rep("female",4),rep("male",1), rep("female",1),rep("male",4), rep("female",2))))
(conditionAdult <- factor(c(rep("cydno", 11), rep("F1", 4))))
(coldata <- data.frame(row.names=colnames(allAdult), conditionAdult,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+conditionAdult)
ddsSex <- DESeq(ddsSex)
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05)
#write.table(resSig, file="Adults_F1vsCYDNO_DE_Cufflinks.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diffAdult<-read.table("Adults_F1vsCYDNO_DE_Cufflinks.txt",header=TRUE)
upcydno<-diffAdult[diffAdult$log2FoldChange >1,]  #
upF1cydno<-diffAdult[diffAdult$log2FoldChange< (-1),] #
all2<-rbind(upcydno,upF1cydno)
all2 <- merge(all2,annotation,by="gene_id")

#SIGNS INVERTED (Deseq puts the group with more data as + fold change)
#CHR18
only18_1_F1<-all2[all2$scaffold == "Hmel218001o",] #
only18_1_F1<-only18_1_F1[only18_1_F1$start >8997,] #
only18_2_F1<-all2[all2$scaffold == "Hmel218002o",] # 
only18_3_F1<-all2[all2$scaffold =="Hmel218003o",] #  
onlypeak_F1<-only18_3_F1[only18_3_F1$end <275070,] # 
onlylodscore_F1<-only18_3_F1[only18_3_F1$end <2393341,] 







##############################       156h APF

###new annotation
#extract from ANNOTATION gene and scaffold info #MELPOMENE
#read annotation #and make tables with gene, scaffold, start, end, INFO
#gff3<-read.delim("combined_156h_Cufflinks_transcripts.gtf", header=F, comment.char="#")
#genesattribute <- data.frame(do.call('rbind', strsplit(as.character(gff3$V9),';',fixed=TRUE)))
#genesattribute2 <- data.frame(do.call('rbind', strsplit(as.character(genesattribute$X1),' ',fixed=TRUE)))
#gff3<-cbind(gff3,genesattribute2)
#annotation<-gff3[c(11,1,4,5)]
#colnames(annotation)[1] <- "gene_id"
#colnames(annotation)[2] <- "scaffold"
#colnames(annotation)[3] <- "start"
#colnames(annotation)[4] <- "end"
#genes<-annotation[!duplicated(annotation$gene_id),]
#write.table(genes, file="156_Cufflinks_gene_info.txt", row.names = FALSE, quote=FALSE)  
annotation<-read.table("156_Cufflinks_gene_info.txt",header=TRUE) 


#read output files from HTseq-count
#mapped to MP
#MPtoMP <-read.table("merged_htseq_156h_APFMPtoCuff_genecounts.txt",header=TRUE)
#CPtoMP <-read.table("merged_htseq_156h_APFCPtoCuff_genecounts.txt",header=TRUE)
#E27toMP <-read.table("merged_htseq_156h_APFE27toCuff_genecounts.txt",header=TRUE)
#all <- merge(MPtoMP,CPtoMP,by="gene_id")
#all <- merge(all,E27toMP,by="gene_id")
#write.table(all, file="156h_toCufflinks_genecounts.txt", row.names = FALSE, quote=FALSE)


#Multi-factor design #sex as random factor   
all156<-read.table("156h_toCufflinks_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(1:19)]
all156<- as.matrix(all156)
#without sex
#(sex <- factor(c(rep("male",1), rep("female", 1),rep("male",2), rep("female",3),rep("male",1), rep("female", 1),rep("male",2), rep("female",4),rep("male",3), rep("female", 1))))   #################
(condition3<- factor(c(rep("mel", 9), rep("cydno", 10))))
(coldata3<- data.frame(row.names=colnames(all156), condition3))
ddsSex3 <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~condition3)  
ddsSex3 <- DESeq(ddsSex3) 
res <- results(ddsSex3)
#resSig <- subset(res, padj < 0.05) 
#write.table(res, file="156h_ALL_toCufflinks.txt", row.names = TRUE, quote=FALSE) 
#write.table(resSig, file="156h_DE_sexfactor_toCufflinks.txt", row.names = TRUE, quote=FALSE) 



diff156<-read.table("156h_DE_sexfactor_toCufflinks.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

#CHR18
only18_1<-all[all$scaffold == "Hmel218001o",] #
only18_1<-only18_1[only18_1$start >8997,] #...
only18_2<-all[all$scaffold == "Hmel218002o",] #
only18_3<-all[all$scaffold =="Hmel218003o",] #  
onlypeak<-only18_3[only18_3$end <275070,] #
onlylodscore<-only18_3[only18_3$end <2393341,] 


###################
#CUFF 13400 = IGluR
# CUFF.13413 = #regucalcin



###################
############## HYBRIDS INTRO LINE
#Multi-factor design #sex as random factor 
all156<-read.table("156h_toCufflinks_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(0,2,3,5,12,14,15,1,4,6,7,8,9,10,11,13,16)] #reorder B LOCUS
#B locus introgressed: 108,116,126,155,183,185 
all156<- as.matrix(all156)
(sex <- factor(c(rep("female",1), rep("male", 6),rep("female",1),rep("male", 2),rep("female",3),rep("male",1),rep("female",1),rep("male",1))))
(condition3 <- factor(c(rep("melchr18", 6), rep("cydno", 10))))
(coldatxa3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156hintrogression_DE_sexfactor_toCuff.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diff156<-read.table("156hintrogression_DE_sexfactor_toCuff.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

#CHR18
only18_1<-all[all$scaffold == "Hmel218001o",] #
only18_1<-only18_1[only18_1$start >8997,] #...
only18_2<-all[all$scaffold == "Hmel218002o",] #
only18_3<-all[all$scaffold =="Hmel218003o",] #  
onlypeak<-only18_3[only18_3$end <275070,] #
onlylodscore<-only18_3[only18_3$end <2393341,] 
#CUFF 134815= IGluR 
#CUFF.13825= regucalcin






##############################       60h APF

###new annotation
#extract from ANNOTATION gene and scaffold info #MELPOMENE
#read annotation #and make tables with gene, scaffold, start, end, INFO
#gff3<-read.delim("combined_60h_Cufflinks_transcripts.gtf", header=F, comment.char="#")
#genesattribute <- data.frame(do.call('rbind', strsplit(as.character(gff3$V9),';',fixed=TRUE)))
#genesattribute2 <- data.frame(do.call('rbind', strsplit(as.character(genesattribute$X1),' ',fixed=TRUE)))
#gff3<-cbind(gff3,genesattribute2)
#annotation<-gff3[c(11,1,4,5)]
#colnames(annotation)[1] <- "gene_id"
#colnames(annotation)[2] <- "scaffold"
#colnames(annotation)[3] <- "start"
#colnames(annotation)[4] <- "end"
#genes<-annotation[!duplicated(annotation$gene_id),]
#write.table(genes, file="60_Cufflinks_gene_info.txt", row.names = FALSE, quote=FALSE)  
annotation<-read.table("60_Cufflinks_gene_info.txt",header=TRUE)


#read output files from HTseq-count
#mapped to MP
#MPtoMP <-read.table("merged_htseq_60h_APF_MPtoCuff_genecounts.txt",header=TRUE)
#CPtoMP <-read.table("merged_htseq_60h_APF_CPtoCuff_genecounts.txt",header=TRUE)
#E27toMP <-read.table("merged_htseq_60h_APF_E27toCuff_genecounts.txt",header=TRUE)
#all <- merge(MPtoMP,CPtoMP,by="gene_id")
#all <- merge(all,E27toMP,by="gene_id")
#write.table(all, file="60h_toCufflinks_genecounts.txt", row.names = FALSE, quote=FALSE)


#Multi-factor design #sex as random factor   
all60<-read.table("60h_toCufflinks_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(1:20)]
all60<- as.matrix(all60)
#(sex <- factor(c(rep("female",1), rep("male", 1),rep("female",1), rep("male",2),rep("female",5),rep("male",2),rep("female",2),rep("male",1),rep("female",2),rep("male",1),rep("female",2))))
(condition60 <- factor(c(rep("mel", 10), rep("cydno", 10))))
#(coldata <- data.frame(row.names=colnames(all60), condition60,sex))
(coldata <- data.frame(row.names=colnames(all60), condition60))
#ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata, design=~sex+condition60)
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata, design=~condition60)
ddsSex <- DESeq(ddsSex)
res <- results(ddsSex)
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.05)
#write.table(res, file="60_nosexfactorALL_toCuff.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
#write.table(resSig, file="60_DE_sexfactor_toCuff.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column


diff60<-read.table("60_DE_sexfactor_toCuff.txt",header=TRUE)
upmel<-diff60[diff60$log2FoldChange >1,]  #
upcydno<-diff60[diff60$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")


#CHR18
only18_1<-all[all$scaffold == "Hmel218001o",] #
only18_1<-only18_1[only18_1$start >8997,] #...
only18_2<-all[all$scaffold == "Hmel218002o",] #
only18_3<-all[all$scaffold =="Hmel218003o",] #  
onlypeak<-only18_3[only18_3$end <275070,] #
onlylodscore<-only18_3[only18_3$end <2393341,] 




######################  E27
#B locus introgressed: 152, 161, 179, 193, 197, 209, 215, 224

#Multi-factor design #sex as random factor 
all60<-read.table("60h_toCufflinks_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(0,1,2,4,8,9,13,16,17,3,5,6,7,10,11,12,14,15)] #reorder
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",1), rep("female", 2),rep("male",1), rep("female",2),rep("male",1), rep("female",5),rep("male",4), rep("female",1))))
(condition60 <- factor(c(rep("melchr18", 8), rep("cydno", 9))))
(coldata <- data.frame(row.names=colnames(all60), condition60,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata, design=~sex+condition60)
ddsSex <- DESeq(ddsSex)
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="60hintrogression_DE_sexfactor_toCuff.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diff60<-read.table("60hintrogression_DE_sexfactor_toCuff.txt",header=TRUE)
upmel<-diff60[diff60$log2FoldChange >1,] #
upcydno<-diff60[diff60$log2FoldChange< (-1),] #
all60<-rbind(upmel,upcydno)
all60 <- merge(all60,annotation,by="gene_id")


#Check at the level of the QTLs    
#CHR18
only18_2<-all60[all60$scaffold == "Hmel218002o",] #
only18_3<-all60[all60$scaffold =="Hmel218003o",] # 3
onlypeak<-only18_3[only18_3$end <275070,] #4
onlylodscore<-only18_3[only18_3$end <2393341,] #


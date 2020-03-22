setwd("/Users/matteo/Desktop/")

library(DESeq2)
library(dplyr)
annotation<-read.table("gene_info.txt",header=TRUE) 

##### CONDUCT DIFFERENTIAL GENE EXPRESSION ANALYSES 

#### between hybrids segregating at the QTL on chromosome 18

#######BC3 hybrids at 156h APF
#B locus introgressed #individuals number 108,116,126,155,183,185 #correspond to column #2,3,7,8,9,11
#Multi-factor design #sex as random factor  
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(0,11,2,3,7,8,9,1,4,5,6,10,12,13,14,15,16)] #reorder
all156<- as.matrix(all156)
(sex <- factor(c(rep("female",1), rep("male", 10),rep("female",5))))
(condition3 <- factor(c(rep("melchr18", 6), rep("cydno", 10))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
#write.table(res, file="156hintrogression_foldchanges_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156hintrogression_DE_sexfactor_toMP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diff156<-read.table("156hintrogression_DE_sexfactor_toMP.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

#Check at the level of the QTLs    
#CHR18
only18_2<-all[all$scaffold == "Hmel218002o",] #
only18_3<-all[all$scaffold =="Hmel218003o",] #
onlylodscore<-only18_3[only18_3$end <2393341,] #

#### USING ONLY MALES
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(0,2,3,7,8,9,1,4,5,6,10)] #reorder
all156<- as.matrix(all156)
(condition3 <- factor(c(rep("melchr18", 5), rep("cydno", 5))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3))
dds <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~condition3)
dds <- DESeq(dds) 
res <- results(dds)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156hintrogression_DE_males_toMP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diff156<-read.table("156hintrogression_DE_males_toMP.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")


#######BC3 hybrids at 60h APF
#B locus introgressed: individuals # 152, 161, 179, 193, 197, 209, 215, 224   #columns 1,2,7,8,10,14,15,17
#Multi-factor design #sex as random factor
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(0,1,2,7,8,10,14,15,17,3,4,5,6,9,11,12,13,16)] #reorder
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",3), rep("female", 5),rep("male",4), rep("female",5))))
(condition60 <- factor(c(rep("melchr18", 8), rep("cydno", 9))))
(coldata <- data.frame(row.names=colnames(all60), condition60,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata, design=~sex+condition60)
ddsSex <- DESeq(ddsSex)
res <- results(ddsSex)
#write.table(res, file="60hintrogression_foldchanges_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="60hintrogression_DE_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diff60<-read.table("60hintrogression_DE_sexfactor.txt",header=TRUE)
upmel<-diff60[diff60$log2FoldChange >1,] #
upcydno<-diff60[diff60$log2FoldChange< (-1),] #
all60<-rbind(upmel,upcydno)
all60 <- merge(all60,annotation,by="gene_id")

#Check at the level of the QTLs    
#CHR18
only18_2<-all60[all60$scaffold == "Hmel218002o",] 
only18_3<-all60[all60$scaffold =="Hmel218003o",] #  
onlylodscore<-only18_3[only18_3$end <2393341,] #

########ONLY MALES
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(0,1,2,7,3,4,5,6)] #reorder
all60<- as.matrix(all60)
(condition60 <- factor(c(rep("melchr18", 3), rep("cydno", 4))))
(coldata <- data.frame(row.names=colnames(all60), condition60))
dds <- DESeqDataSetFromMatrix(countData=all60, colData=coldata, design=~condition60)
dds <- DESeq(dds)
res <- results(dds)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="60hintrogression_DE_males_toMP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diff60<-read.table("60hintrogression_DE_males_toMP.txt",header=TRUE)
upmel<-diff60[diff60$log2FoldChange >1,] #
upcydno<-diff60[diff60$log2FoldChange< (-1),] #
all60<-rbind(upmel,upcydno)
all60 <- merge(all60,annotation,by="gene_id")



######  CHECK DIFFERENTIAL GENE EXPRESSION WHEN COMPARING BC3 HYBRIDS SEGREGATING AT OTHER CHROMOSOMES

################## CHR 4
###156 h  APF   
#### chr 4 introgressed: 108, 116, 126, 136, 139, 154, 183, 185, 189
#### chr 4 not introgressed: 105, 123, 131, 133, 140, 155, 166
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(0,2,3,6,8,9,10,11,13,14,1,4,5,7,12,15,16)] #reorder
all156<- as.matrix(all156)
(sex <- factor(c(rep("male",6), rep("female",3 ),rep("male",4), rep("female",3 ))))
(condition3 <- factor(c(rep("melchr4", 9), rep("cydno", 7))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156hintrogression_CHR4_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("156hintrogression_CHR4_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

#### 60 h APF
#### chr 4 introgressed: 152, 165, 179, 187, 188, 193, 197, 198, 199, 212, 214, 224
#### chr 4 not introgressed: 161, 192, 201, 209, 215
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(0,1,2,3,4,6,9,10,11,12,14,16,17,5,7,8,13,15)] #reorder 
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",5), rep("female",7 ),rep("male",2), rep("female",3 ))))
(condition3 <- factor(c(rep("melchr4", 12), rep("cydno", 5))))
(coldata3 <- data.frame(row.names=colnames(all60), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="60hintrogression_CHR4_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("60hintrogression_CHR4_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")


###### CHR 1
#### 156 h APF
#### chr 1 introgressed: 105, 126, 136, 139, 140, 154, 155, 183, 185, 189
#### chr 1 not introgressed: 108, 116, 123, 131, 133, 166
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(0,1,3,6,7,8,9,10,13,14,15,2,4,5,11,12,16)] #reorder
all156<- as.matrix(all156)
(sex <- factor(c(rep("male",7), rep("female",3 ),rep("male",3), rep("female",3 ))))
(condition3 <- factor(c(rep("melchr1", 10), rep("cydno", 6))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156hintrogression__CHR1_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("156hintrogression__CHR1_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

#### 60 h APF
#### chr 1 introgressed: 152, 161, 193, 198, 201, 212, 224
#### chr 1 not introgressed:  165, 179, 187, 188, 192, 197, 199, 209, 214, 215
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(0,1,2,3,5,12,13,15,4,6,7,8,9,10,11,14,16,17 )] #reorder #first chr1 introgressed then no
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",4), rep("female",1 ),rep("male",1), rep("female",1 ),rep("male",1), rep("female",6 ),rep("male",1), rep("female",2 ))))
(condition3 <- factor(c(rep("melchr1", 7), rep("cydno", 10))))
(coldata3 <- data.frame(row.names=colnames(all60), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="60hintrogression_CHR1_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("60hintrogression_CHR1_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")


#####  CHR 15
###156 h APF  
#### chr 15 introgressed: 116, 123, 126, 131, 136, 154, 155, 183, 185, 189
#### chr 15 not introgressed: 105, 108, 133, 139, 140, 166
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(2,3,4,6,7,8,9,10,12,13,1,5,11,14,15,16)] #reorder
all156<- as.matrix(all156)
(sex <- factor(c(rep("male",8), rep("female",2 ),rep("male",2), rep("female",4))))
(condition3 <- factor(c(rep("melchr15", 10), rep("cydno", 6))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156hintrogression__CHR15_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("156hintrogression__CHR15_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")
###60 h APF     
#### chr 15 introgressed: 161, 179, 192, 193, 197, 199, 209, 214
#### chr 15 not introgressed: 152, 165, 187, 188, 198, 201, 212, 215, 224
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(2,4,8,10,13,14,15,16,1,3,5,6,7,9,11,12,17)] #reorder
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",2), rep("female",6),rep("male",5), rep("female",4))))
(condition3 <- factor(c(rep("melchr15", 8), rep("cydno", 9))))
(coldata3 <- data.frame(row.names=colnames(all60), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="60hintrogression_CHR15_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("60hintrogression_CHR15_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #3
upcydno<-diff156[diff156$log2FoldChange< (-1),] #1
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")


#### CHR 20
###156 h   APF  
#### chr 20 introgressed: 105, 123, 131, 133, 140, 154, 166, 183, 189
#### chr 20 not introgressed: 108, 116, 126, 136, 139, 155, 185
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(1,4,5,6,8,10,12,15,16,2,3,7,9,11,13,14)] #reorder
all156<- as.matrix(all156)
(sex <- factor(c(rep("male",6), rep("female",3 ),rep("male",4), rep("female",3))))
(condition3 <- factor(c(rep("melchr20", 9), rep("cydno", 7))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
#write.table(res, file="156hintrogression_CHR20ALL_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156hintrogression__CHR20_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("156hintrogression__CHR20_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

###60 h APF
#### chr 20 introgressed: 152, 161, 187, 188, 192, 193, 198, 209, 212, 214, 215, 224
#### chr 20 not introgressed: 165, 179, 197, 199, 201
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(1,2,3,6,7,8,11,12,13,15,16,17,4,5,9,10,14)] #reorder
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",5), rep("female",7),rep("male",2), rep("female",3))))
(condition3 <- factor(c(rep("melchr15", 12), rep("cydno", 5))))
(coldata3 <- data.frame(row.names=colnames(all60), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="60hintrogression_CHR20_sexfactor.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("60hintrogression_CHR20_sexfactor.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

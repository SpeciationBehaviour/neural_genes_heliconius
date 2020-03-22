setwd("/Users/matteo/Desktop/")  

library(DESeq2)
library(dplyr)

#CONDUCT DIFFERENTIAL GENE EXPRESSION ANALYSES USING THE TRRANSCRIPT-BASED ANNOTATION

###Example with Adult species samples/comparisons

#extract from ANNOTATION gene and scaffold info 
#read annotation #and make tables with gene, scaffold, start, end, INFO
gff3<-read.delim("Adults_combined_Cufflinks_transcripts.gtf", header=F, comment.char="#")
genesattribute <- data.frame(do.call('rbind', strsplit(as.character(gff3$V9),';',fixed=TRUE)))
genesattribute2 <- data.frame(do.call('rbind', strsplit(as.character(genesattribute$X1),' ',fixed=TRUE)))
gff3<-cbind(gff3,genesattribute2)
annotation<-gff3[c(11,1,4,5)]
colnames(annotation)[1] <- "gene_id"
colnames(annotation)[2] <- "scaffold"
colnames(annotation)[3] <- "start"
colnames(annotation)[4] <- "end"
genes<-annotation[!duplicated(annotation$gene_id),]
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

# Melpomene vs Cydno
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
#write.table(resSig, file="Adults_DE_sexfactor_toCufflinks.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
#write.table(res, file="Adults_sexfactorALL_toCufflinks.txt", row.names = TRUE, quote=FALSE) 

#Differentially expressed genes:
diffAdult<-read.table("Adults_DE_sexfactor_toCufflinks.txt",header=TRUE)
upmel<-diffAdult[diffAdult$log2FoldChange >1,]  #
upcydno<-diffAdult[diffAdult$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

#See which differentially expresse genes at the level of the QTL on CHR18
only18_1<-all[all$scaffold == "Hmel218001o",] 
only18_1<-only18_1[only18_1$start >8997,] #...
only18_2<-all[all$scaffold == "Hmel218002o",] 
only18_3<-all[all$scaffold =="Hmel218003o",] #  
onlypeak<-only18_3[only18_3$end <275070,] #  
onlylodscore<-only18_3[only18_3$end <2393341,] 
#CUFF.13201 #Ionotropic glut receptor
##CUFF.13212 #regucalcin2


###########     F1 vs cydno
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

#Check which genes at the level of the QTL on CHR18
#NOTE: SIGNS INVERTED (DESeq2 > the group with more data as + fold change)
only18_1_F1<-all2[all2$scaffold == "Hmel218001o",] #
only18_1_F1<-only18_1_F1[only18_1_F1$start >8997,] #
only18_2_F1<-all2[all2$scaffold == "Hmel218002o",] # 
only18_3_F1<-all2[all2$scaffold =="Hmel218003o",] #  
onlypeak_F1<-only18_3_F1[only18_3_F1$end <275070,] # 
onlylodscore_F1<-only18_3_F1[only18_3_F1$end <2393341,] 

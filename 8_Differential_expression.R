#Conduct differential gene expression analyses between H. melpomene and H.cydno, across developmental stages

setwd("/Users/matteo/Desktop/")

#install DESeq2
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
library(dplyr)

#extract from Hmel2.5 ANNOTATION gene and scaffold info 
#download annotation at:http://butterflygenome.org/?q=node/4
#read annotation #and make tables with gene, scaffold, start, end, INFO
gff3<-read.delim("Hmel2.5.gff3", header=F, comment.char="#")
gff.genes <- gff3[gff3[,3]=="gene",]
genesattribute <- data.frame(do.call('rbind', strsplit(as.character(gff.genes$V9),';',fixed=TRUE)))
genesattribute2 <- data.frame(do.call('rbind', strsplit(as.character(genesattribute$X1),'=',fixed=TRUE)))
gff.genes2<-cbind(gff.genes,genesattribute2)
annotation<-gff.genes2[c(11,1,4,5)]
colnames(annotation)[1] <- "gene_id"
colnames(annotation)[2] <- "scaffold"
colnames(annotation)[3] <- "start"
colnames(annotation)[4] <- "end"
#write.table(annotation, file="gene_info.txt", row.names = FALSE, quote=FALSE)  
annotation<-read.table("gene_info.txt",header=TRUE) 

#extract from ANNOTATION gene and scaffold info #CYDNO
#read annotation #and make tables with gene, scaffold, start, end, INFO
#gff3<-read.delim("Hcyd_corrected.gff3", header=F, comment.char="#")
#gff.genes <- gff3[gff3[,3]=="gene",]
#genesattribute <- data.frame(do.call('rbind', strsplit(as.character(gff.genes$V9),';',fixed=TRUE)))
#genesattribute2 <- data.frame(do.call('rbind', strsplit(as.character(genesattribute$X1),'=',fixed=TRUE)))
#gff.genes2<-cbind(gff.genes,genesattribute2)
#annotation<-gff.genes2[c(11,1,4,5)]
#colnames(annotation)[1] <- "gene_id"
#colnames(annotation)[2] <- "scaffold"
#colnames(annotation)[3] <- "start"
#colnames(annotation)[4] <- "end"
#write.table(annotation, file="gene_info_cydno.txt", row.names = FALSE, quote=FALSE)  
annotation_cyd<-read.table("gene_info_cydno.txt",header=TRUE)

#Just looking at genes at the level of the QTLs
#CHR 18
#PEAK #min 8997 to end - Hmel218001o  #whole Hmel218002o #from 0 to 275070 - Hmel218003o
#LOD 1.5 score same as above but until 4409863 Hmel218003o
only18_1<-annotation[annotation$scaffold == "Hmel218001o",] 
only18_1<-only18_1[only18_1$start >8997,] #
only18_2<-annotation[annotation$scaffold == "Hmel218002o",] #
only18_3<-annotation[annotation$scaffold =="Hmel218003o",]  
onlylodscore18<-only18_3[only18_3$end <2393341,] #149
#200 genes in total

#READ ADULTS SAMPLES
#read output files from HTseq-count AND merge in one table
MPtoMP <-read.table("merged_htseq_MPtoMP_genecounts.txt",header=TRUE)
CPtoMP <-read.table("merged_htseq_CPtoMP_genecounts.txt",header=TRUE)
F1toMP <-read.table("merged_htseq_F1toMP_genecounts.txt",header=TRUE)
all <- merge(MPtoMP,CPtoMP,by="gene_id")
all <- merge(all,F1toMP,by="gene_id")
#write.table(all, file="Adults_toMP_genecounts.txt", row.names = FALSE, quote=FALSE)
#(do the same for samples at 156h after pupal formation (APF) and 60h APF)


######### CONDUCT DIFFERENTIAL GENE EXPRESSION ANALYSES with DESeq2

######### Adults comparisons
#Multi-factor design #sex as random factor   
allAdult<-read.table("Adults_toMP_gene_counts.txt",header=TRUE, row.names = 1)
allAdult<-allAdult[,c(1:23)]  
allAdult<- as.matrix(allAdult)
(sex <- factor(c(rep("male",7), rep("female", 5),rep("male",4), rep("female",7))))
(conditionAdult <- factor(c(rep("mel", 12), rep("cydno", 11))))
(coldata <- data.frame(row.names=colnames(allAdult), conditionAdult,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+conditionAdult)
ddsSex <- DESeq(ddsSex)
res <- results(ddsSex)
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.05) #restrict to significant 
#write.table(resSig, file="Adults_DE_sexfactor_toMP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
#write.table(res, file="Adults_sexfactor_all.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

#Considering differentially expressed those showing a 2fold-change in expression besides padj < 0.05
diffAdult<-read.table("Adults_DE_sexfactor_toMP.txt",header=TRUE)
upmel<-diffAdult[diffAdult$log2FoldChange >1,]  
upcydno<-diffAdult[diffAdult$log2FoldChange< (-1),] 
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

#Now check at the level of the QTL on CHR 18
only18_1<-all[all$scaffold == "Hmel218001o",] 
only18_1<-only18_1[only18_1$start >8997,] 
only18_2<-all[all$scaffold == "Hmel218002o",]  #9
only18_3<-all[all$scaffold =="Hmel218003o",]   
onlypeak<-only18_3[only18_3$end <275070,] #   
onlylodscore<-only18_3[only18_3$end <2393341,] #11


#Check which genes are differentially expressed when mapping to the H.cydno genome (above is with melpomene)
#Multi-factor design #sex as random factor 
allAdult2<-read.table("Adults_toCP_gene_counts.txt",header=TRUE, row.names = 1)
allAdult2<-allAdult2[,c(1:23)]  
allAdult2<- as.matrix(allAdult2)
(sex <- factor(c(rep("male",7), rep("female", 5),rep("male",4), rep("female",7))))
(conditionAdult2 <- factor(c(rep("mel", 12), rep("cydno", 11))))
(coldata2 <- data.frame(row.names=colnames(allAdult2), conditionAdult2,sex))
ddsSex2 <- DESeqDataSetFromMatrix(countData=allAdult2, colData=coldata2, design=~sex+conditionAdult2)
dds2 <- DESeq(ddsSex2)
res <- results(dds2)
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.05)
#write.table(resSig, file="Adults_DE_sexfactor_toCP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diffAdult_cyd<-read.table("Adults_DE_sexfactor_toCP.txt",header=TRUE)
upmel_cyd<-diffAdult_cyd[diffAdult_cyd$log2FoldChange >1,]  #
upcydno_cyd<-diffAdult_cyd[diffAdult_cyd$log2FoldChange< (-1),] #
all_cyd<-rbind(upmel_cyd,upcydno_cyd)
all_cyd <- merge(all_cyd,annotation_cyd,by="gene_id")

#SEE OVERLAP GENES CONSIDERED DIFFERENTIALLY EXPRESSED when mapping to both genomes at the QTL on CHR 18
overlap_peak_scf2<-inner_join(only18_2,all_cyd,by="gene_id") #
overlap_lodscore<-inner_join(onlylodscore,all_cyd,by="gene_id") #


######## Pupae 156h APF comparisons
#Multi-factor design #sex as random factor   
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(1:19)]
all156<- as.matrix(all156)
(sex <- factor(c(rep("male",4), rep("female", 5),rep("male",5), rep("female",5))))
(condition3 <- factor(c(rep("mel", 9), rep("cydno", 10))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156h_DE_sexfactor_toMP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
#write.table(res, file="156h_sexfactor_all.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column
diff156<-read.table("156h_DE_sexfactor_toMP.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")
#QTL on CHR 18
only18_1<-all[all$scaffold == "Hmel218001o",] #
only18_1<-only18_1[only18_1$start >8997,] #...
only18_2<-all[all$scaffold == "Hmel218002o",] #
only18_3<-all[all$scaffold =="Hmel218003o",] #   
onlypeak<-only18_3[only18_3$end <275070,]  
onlylodscore<-only18_3[only18_3$end <2393341,] #

####### MAPPING TO CYDNO
##Multi-factor design #sex as random factor   
all156<-read.table("156h_toCP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(1:19)]
all156<- as.matrix(all156)
(sex <- factor(c(rep("male",4), rep("female", 5),rep("male",5), rep("female",5))))
(condition3 <- factor(c(rep("mel", 9), rep("cydno", 10))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex <- DESeq(ddsSex) 
res <- results(ddsSex)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="156h_DE_sexfactor_toCP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diff156CP<-read.table("156h_DE_sexfactor_toCP.txt",header=TRUE)
upmelCP<-diff156CP[diff156CP$log2FoldChange >1,] #518
upcydnoCP<-diff156CP[diff156CP$log2FoldChange< (-1),] #403
allCP<-rbind(upmelCP,upcydnoCP)
all_CP<- merge(allCP,annotation_cyd,by="gene_id") #921

#OVERLAP GENES DIFF EXPRESSED at QTL chr 18
overlap_peak_scf2<-inner_join(only18_2,all_CP,by="gene_id") #
overlap_lodscore<-inner_join(onlylodscore,all_CP,by="gene_id") #



## See if the mapping genome introduces biases in differential expression inference (e.g. when mapping to the melpomene genome, more melpomene up-regulated genes? artifact of melpomene RNA reads that mapped better?)
#FISHER-TEST
#ADULTS
up_regulated <-
  matrix(c(694, 733, 390, 451),
         nrow = 2,
         dimnames = list(mel = c("up_mel", "up_cyd"),
                         cydn = c("up_mel", "up_cyd")))
fisher.test(up_regulated) #two sided

#156hAPF
up_regulated <-
  matrix(c(837, 667, 518, 403),
         nrow = 2,
         dimnames = list(mel = c("up_mel", "up_cyd"),
                         cydn = c("up_mel", "up_cyd")))
fisher.test(up_regulated)

#60hAPF
up_regulated <-
  matrix(c(846, 642, 490, 376),
         nrow = 2,
         dimnames = list(mel = c("up_mel", "up_cyd"),
                         cydn = c("up_mel", "up_cyd")))
fisher.test(up_regulated) 



###See if chromosome 18 is enriched for differentially expressed genes compared to other chromosome
#Fisher's exact test. 
#(1) The number of differentially expressed genes on chromosome 18, 
#(2) the number of non-DE genes on chromosome 18 
#(3) the number of DE genes not on chromosome other autosomes
#(4) the number of non-DE genes not on other autosomes

#Adults
annotation<-read.table("chromo_gene_info.txt",header=TRUE)
diffAdult<-read.table("Adults_DE_sexfactor_toMP.txt",header=TRUE)
upmel<-diffAdult[diffAdult$log2FoldChange >1,]  
upcydno<-diffAdult[diffAdult$log2FoldChange< (-1),] 
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

chr18_DE<-all[all$scaffold=="chr18",] #75
chr18<-annotation[annotation$scaffold=="chr18",] #992
chr18_notDE<-anti_join(chr18,chr18_DE,by="gene_id") #917

otherautosomes_DE<-all[!(all$scaffold=="chr18") & !(all$scaffold=="chr21"),] #1199
otherautosomes<-annotation[!(annotation$scaffold=="chr18") & !(annotation$scaffold=="chr21"),] #17719
otherautosomes_notDE<-anti_join(otherautosomes,otherautosomes_DE,by="gene_id") #16520

DE <-
  matrix(c(75, 917, 1199, 16520),
         nrow = 2,
         dimnames = list(chr18 = c("up_18", "not_18"),
                         autosomes = c("up_autosm", "not_autosm")))
fisher.test(DE) #two sided
#p-value = 0.3312


#156h APF
diff156<-read.table("156h_DE_sexfactor_toMP.txt",header=TRUE)
upmel<-diff156[diff156$log2FoldChange >1,] #
upcydno<-diff156[diff156$log2FoldChange< (-1),] #
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

chr18_DE<-all[all$scaffold=="chr18",] #61
chr18_notDE<-anti_join(chr18,chr18_DE,by="gene_id") #931

otherautosomes_DE<-all[!(all$scaffold=="chr18") & !(all$scaffold=="chr21"),] #1277
otherautosomes_notDE<-anti_join(otherautosomes,otherautosomes_DE,by="gene_id") #16442

DE <-
  matrix(c(61, 931, 1277, 16442),
         nrow = 2,
         dimnames = list(chr18 = c("up_18", "not_18"),
                         autosomes = c("up_autosm", "not_autosm")))
fisher.test(DE) #two sided
#p-value = 0.2288



#60 h APF
diff60<-read.table("60_DE_sexfactor_toMP.txt",header=TRUE)
upmel<-diff60[diff60$log2FoldChange >1,]  
upcydno<-diff60[diff60$log2FoldChange< (-1),] 
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")

chr18_DE<-all[all$scaffold=="chr18",] #76
chr18_notDE<-anti_join(chr18,chr18_DE,by="gene_id") #916

otherautosomes_DE<-all[!(all$scaffold=="chr18") & !(all$scaffold=="chr21"),] #1262
otherautosomes_notDE<-anti_join(otherautosomes,otherautosomes_DE,by="gene_id") #16457

DE <-
  matrix(c(76, 916, 1262, 16457),
         nrow = 2,
         dimnames = list(chr18 = c("up_18", "not_18"),
                         autosomes = c("up_autosm", "not_autosm")))
fisher.test(DE) #two sided
#p-value = 0.5264


#QTl on chr 18 , rest of the genes
#all - QTl18 = rest
#Adults: 1374 - 20 = 1354
#156h: 1504 -14= 1490
#60h: 1488 -15 = 1473

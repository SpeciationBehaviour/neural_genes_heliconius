#Note that, as in	S. H. Martin, J. W. Davey, C. Salazar, C. D. Jiggins, Recombination rate variation shapes barriers to introgression across butterfly genomes. PLoS Biol. 17, 1–28 (2019), the genome assembly used is in Hmel2 scaffolds, but the coordinate system has been converted to chromosomes according to the Hmel2.5 scaffolding
#Two distinct graphs for scaffold 2 and scaffold 3 of chromosome 18 are produced and then merged 
#Note that there is no evidence for chromosomal inversions (longer than 50kb) in melpomene and cydno (Davey, J. W., Barker, S. L., Rastas, P. M., Pinharanda, A., Martin, S. H., Durbin, R., ... & Jiggins, C. D. (2017). No evidence for maintenance of a sympatric Heliconius species barrier by chromosomal inversions. Evolution Letters, 1(3), 138-154.). The gap between scaffold 2 and 3 is likely to contain repeat elements only, and should not be indicative of an inversion.

setwd("/Users/matteo/Desktop/")

#source("https://bioconductor.org/biocLite.R")
#biocLite("Gviz")
#biocLite("DESeq2")
#install.packages('stringi')
library("Gviz")
library("DESeq2")
options(ucscChromosomeNames=FALSE)

#make genomic ranges #to later match genes physical position and fold-change in expression level
annotation<-read.table("gene_info.txt",header=TRUE) 
#in dds2 results ordered by gene name, so order gene ranges also by gene name
annotation <- annotation[order(annotation$gene_id),]
GR<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))
a <- AnnotationTrack(GR, name = "gene ranges") #create object that stores genomic coordinates for every gene

### RECONDUCT DIFFERENTIAL GENE EXPRESSION ANALYSES and store values
#SPECIES
### Adult stage
allAdult<-read.table("Adults_toMP_gene_counts.txt",header=TRUE, row.names = 1)
allAdult<-allAdult[,c(1:23)]  
allAdult<- as.matrix(allAdult)
(sex <- factor(c(rep("male",7), rep("female", 5),rep("male",4), rep("female",7))))
(conditionAdult <- factor(c(rep("mel", 12), rep("cydno", 11))))
(coldata <- data.frame(row.names=colnames(allAdult), conditionAdult,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+conditionAdult)
rowRanges(ddsSex) <- GR #add genomic ranges to differential expression results
dds2 <- DESeq(ddsSex)
resGR <- results(dds2, format = "GRanges")
resGR$padj<--log10(resGR$padj) #transform p values adjusted in -log base10
d <- DataTrack(resGR, data = "padj", type = "p", name = "-log10(padj)", strand = "*") #store p-values and genomic ranges in Gvizobject

#### 156hAPF
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(1:19)]
all156<- as.matrix(all156)
(sex <- factor(c(rep("male",4), rep("female", 5),rep("male",5), rep("female",5))))
(condition3 <- factor(c(rep("mel", 9), rep("cydno", 10))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex2 <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
ddsSex2 <- DESeq(ddsSex2) 
rowRanges(ddsSex2) <- GR
dds2 <- DESeq(ddsSex2)
resGR3 <- results(dds2, format = "GRanges")
resGR3$padj<--log10(resGR3$padj) 
d3 <- DataTrack(resGR3, data = "padj", type = "p", name = "-log10(padj)", strand = "*")

### 60hAPF
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(1:20)]
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",3), rep("female", 7),rep("male",4), rep("female",6))))
(condition60 <- factor(c(rep("mel", 10), rep("cydno", 10))))
(coldata <- data.frame(row.names=colnames(all60), condition60,sex))
ddsSex3 <- DESeqDataSetFromMatrix(countData=all60, colData=coldata, design=~sex+condition60)
ddsSex3<- DESeq(ddsSex3)
rowRanges(ddsSex3) <- GR
dds2 <- DESeq(ddsSex3)
resGR4 <- results(dds2, format = "GRanges")
resGR4$padj<--log10(resGR4$padj) #log10(0.05) #1.3
d4 <- DataTrack(resGR4, data = "padj", type = "p", name = "-log10(padj)", strand = "*")

## store values with associated graphical parameters
d <- DataTrack(resGR, data = "padj", type = "p", name = "-log10(padj)", strand = "*",background.title = "white",col.axis=1,col.title=1,col="black",baseline=1.3, col.baseline="black",ylim=c(0,30))
d3 <- DataTrack(resGR3, data = "padj", type = "p", name = "-log10(padj)", strand = "*",background.title = "white",col.axis=1,col.title=1,col="black",baseline=1.3, col.baseline="black",ylim=c(0,30))
d4 <- DataTrack(resGR4, data = "padj", type = "p", name = "-log10(padj)", strand = "*",background.title = "white",col.axis=1,col.title=1,col="black",baseline=1.3, col.baseline="black",ylim=c(0,30))

#Set active Chromosome and then withdraw p values
chromosome(d) <- 'Hmel218002o' 
chromosome(d3) <- 'Hmel218002o' 
chromosome(d4) <- 'Hmel218002o' 
p_raw_d_Hmel218002o<-values(d)
p_raw_d3_Hmel218002o<-values(d3)
p_raw_d4_Hmel218002o<-values(d4)

chromosome(d) <- 'Hmel218003o' 
chromosome(d3) <- 'Hmel218003o' 
chromosome(d4) <- 'Hmel218003o' 
p_raw_d_Hmel218003o<-values(d)
p_raw_d3_Hmel218003o<-values(d3)
p_raw_d4_Hmel218003o<-values(d4)

#Set safe range
write.table(d@range,"d.txt")

#Read in ranges and change all factor columns to character columns.
d_seq<-read.table("d.txt",header=T)
for(i in 1:length(d_seq[1,])){
  if(is.factor(d_seq[,i])){
    d_seq[,i]<-levels(d_seq[,i])[d_seq[,i]]
  }
}

#Repeat the same for other two
write.table(d3@range,"d3.txt")

d3_seq<-read.table("d.txt",header=T)
for(i in 1:length(d3_seq[1,])){
  if(is.factor(d3_seq[,i])){
    d3_seq[,i]<-levels(d3_seq[,i])[d3_seq[,i]]
  }
}

write.table(d4@range,"d4.txt")

d4_seq<-read.table("d.txt",header=T)
for(i in 1:length(d4_seq[1,])){
  if(is.factor(d4_seq[,i])){
    d4_seq[,i]<-levels(d4_seq[,i])[d4_seq[,i]]
  }
}

#Retrieve chromosome names
seqnames_d<-d_seq$seqnames
seqnames_d3<-d3_seq$seqnames
seqnames_d4<-d4_seq$seqnames

#Retrieve mean position of scaffold (mean of beginning and end of scaffold)
mean_pos_d<-(d_seq$start+d_seq$end)/2
mean_pos_d3<-(d3_seq$start+d4_seq$end)/2
mean_pos_d4<-(d4_seq$start+d3_seq$end)/2

#Extract those positions that belong to the two scaffolds you're interested in.
mean_pos_d_Hmel218002o<-mean_pos_d[seqnames_d=="Hmel218002o"]
mean_pos_d3_Hmel218002o<-mean_pos_d3[seqnames_d3=="Hmel218002o"]
mean_pos_d4_Hmel218002o<-mean_pos_d4[seqnames_d4=="Hmel218002o"]

mean_pos_d_Hmel218003o<-mean_pos_d[seqnames_d=="Hmel218003o"]
mean_pos_d3_Hmel218003o<-mean_pos_d3[seqnames_d3=="Hmel218003o"]
mean_pos_d4_Hmel218003o<-mean_pos_d4[seqnames_d4=="Hmel218003o"]

#Fold changes
fold_d_Hmel218002o<-resGR$log2FoldChange[seqnames_d=="Hmel218002o"]
fold_d3_Hmel218002o<-resGR3$log2FoldChange[seqnames_d=="Hmel218002o"]
fold_d4_Hmel218002o<-resGR4$log2FoldChange[seqnames_d=="Hmel218002o"]

fold_d_Hmel218003o<-resGR$log2FoldChange[seqnames_d=="Hmel218003o"]
fold_d3_Hmel218003o<-resGR3$log2FoldChange[seqnames_d=="Hmel218003o"]
fold_d4_Hmel218003o<-resGR4$log2FoldChange[seqnames_d=="Hmel218003o"]

#Which positions are free of NA
no_NA_d_Hmel218002o<-sapply(1:length(fold_d_Hmel218002o),function(x) sum(is.na(p_raw_d_Hmel218002o[x]),is.na(mean_pos_d_Hmel218002o[x]),is.na(fold_d_Hmel218002o[x])))
no_NA_d3_Hmel218002o<-sapply(1:length(fold_d3_Hmel218002o),function(x) sum(is.na(p_raw_d3_Hmel218002o[x]),is.na(mean_pos_d3_Hmel218002o[x]),is.na(fold_d3_Hmel218002o[x])))
no_NA_d4_Hmel218002o<-sapply(1:length(fold_d4_Hmel218002o),function(x) sum(is.na(p_raw_d4_Hmel218002o[x]),is.na(mean_pos_d4_Hmel218002o[x]),is.na(fold_d4_Hmel218002o[x])))

no_NA_d_Hmel218003o<-sapply(1:length(fold_d_Hmel218003o),function(x) sum(is.na(p_raw_d_Hmel218003o[x]),is.na(mean_pos_d_Hmel218003o[x]),is.na(fold_d_Hmel218003o[x])))
no_NA_d3_Hmel218003o<-sapply(1:length(fold_d3_Hmel218003o),function(x) sum(is.na(p_raw_d3_Hmel218003o[x]),is.na(mean_pos_d3_Hmel218003o[x]),is.na(fold_d3_Hmel218003o[x])))
no_NA_d4_Hmel218003o<-sapply(1:length(fold_d4_Hmel218003o),function(x) sum(is.na(p_raw_d4_Hmel218003o[x]),is.na(mean_pos_d4_Hmel218003o[x]),is.na(fold_d4_Hmel218003o[x])))

#Color vectors
col_d_Hmel218002o_alt<-sapply((1:length(fold_d_Hmel218002o))[no_NA_d_Hmel218002o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d_Hmel218002o[x]<=1.3|(fold_d_Hmel218002o[x]>=(-1)&fold_d_Hmel218002o[x]<=1),p_raw_d_Hmel218002o[x]>1.3&fold_d_Hmel218002o[x]<(-1),p_raw_d_Hmel218002o[x]>1.3&fold_d_Hmel218002o[x]>1)])
col_d3_Hmel218002o_alt<-sapply((1:length(fold_d3_Hmel218002o))[no_NA_d3_Hmel218002o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d3_Hmel218002o[x]<=1.3|(fold_d3_Hmel218002o[x]>=(-1)&fold_d3_Hmel218002o[x]<=1),p_raw_d3_Hmel218002o[x]>1.3&fold_d3_Hmel218002o[x]<(-1),p_raw_d3_Hmel218002o[x]>1.3&fold_d3_Hmel218002o[x]>1)])
col_d4_Hmel218002o_alt<-sapply((1:length(fold_d4_Hmel218002o))[no_NA_d4_Hmel218002o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d4_Hmel218002o[x]<=1.3|(fold_d4_Hmel218002o[x]>=(-1)&fold_d4_Hmel218002o[x]<=1),p_raw_d4_Hmel218002o[x]>1.3&fold_d4_Hmel218002o[x]<(-1),p_raw_d4_Hmel218002o[x]>1.3&fold_d4_Hmel218002o[x]>1)])

col_d_Hmel218003o_alt<-sapply((1:length(fold_d_Hmel218003o))[no_NA_d_Hmel218003o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d_Hmel218003o[x]<=1.3|(fold_d_Hmel218003o[x]>=(-1)&fold_d_Hmel218003o[x]<=1),p_raw_d_Hmel218003o[x]>1.3&fold_d_Hmel218003o[x]<(-1),p_raw_d_Hmel218003o[x]>1.3&fold_d_Hmel218003o[x]>1)])
col_d3_Hmel218003o_alt<-sapply((1:length(fold_d3_Hmel218003o))[no_NA_d3_Hmel218003o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d3_Hmel218003o[x]<=1.3|(fold_d3_Hmel218003o[x]>=(-1)&fold_d3_Hmel218003o[x]<=1),p_raw_d3_Hmel218003o[x]>1.3&fold_d3_Hmel218003o[x]<(-1),p_raw_d3_Hmel218003o[x]>1.3&fold_d3_Hmel218003o[x]>1)])
col_d4_Hmel218003o_alt<-sapply((1:length(fold_d4_Hmel218003o))[no_NA_d4_Hmel218003o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d4_Hmel218003o[x]<=1.3|(fold_d4_Hmel218003o[x]>=(-1)&fold_d4_Hmel218003o[x]<=1),p_raw_d4_Hmel218003o[x]>1.3&fold_d4_Hmel218003o[x]<(-1),p_raw_d4_Hmel218003o[x]>1.3&fold_d4_Hmel218003o[x]>1)])


### HYBRIDS COMPARISONS
#F1 VS CYDNO
allAdult<-read.table("Adults_toMP_gene_counts.txt",header=TRUE, row.names = 1)
allAdult<-allAdult[,c(13:24,26:28)] 
allAdult<-allAdult[,c(12:15,0:11)] 
allAdult<- as.matrix(allAdult)
(sex <- factor(c(rep("male",2), rep("female", 2),rep("male",4), rep("female",7))))
(conditionAdult <- factor(c(rep("F1", 4), rep("cydno", 11))))
(coldata <- data.frame(row.names=colnames(allAdult), conditionAdult,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+conditionAdult)
rowRanges(ddsSex ) <- GR
ddsSex  <- DESeq(ddsSex)
results<-results(ddsSex) 
results #fold changes reversed #now negative means more in F1 (change later)
resGR21 <- results(ddsSex , format = "GRanges")
resGR21$padj<--log10(resGR21$padj)
d21 <- DataTrack(resGR21, data = "padj", type = "p", name = "-log10(padj)", strand = "*",col=c("black"),background.title = "white",col.axis=1,col.title=1,ylim=c(0,25))

#INTROGRESSION LINE 156h APF
all156<-read.table("156h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all156<-all156[,c(0,20:35)]
all156<-all156[,c(0,11,2,3,7,8,9,1,4,5,6,10,12,13,14,15,16)] #reorder
all156<- as.matrix(all156)
(sex <- factor(c(rep("female",1), rep("male", 10),rep("female",5))))
(condition3 <- factor(c(rep("melchr18", 6), rep("cydno", 10))))
(coldata3 <- data.frame(row.names=colnames(all156), condition3,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all156, colData=coldata3, design=~sex+condition3)
rowRanges(ddsSex) <- GR
ddsSex <- DESeq(ddsSex) 
resGR9 <- results(ddsSex, format = "GRanges")
resGR9$padj<--log10(resGR9$padj)
d9 <- DataTrack(resGR9, data = "padj", type = "p", name = "-log10(padj)", strand = "*",col=c("black"),background.title = "white",col.axis=1,col.title=1,ylim=c(0,25))

#INTROGRESSION LINE 60h APF
all60<-read.table("60h_toMP_genecounts.txt",header=TRUE, row.names = 1)
all60<-all60[,c(0,21:37)]
all60<-all60[,c(0,1,2,7,8,10,14,15,17,3,4,5,6,9,11,12,13,16)] #reorder
all60<- as.matrix(all60)
(sex <- factor(c(rep("male",3), rep("female", 5),rep("male",4), rep("female",5))))
(condition60 <- factor(c(rep("melchr18", 8), rep("cydno", 9))))
(coldata <- data.frame(row.names=colnames(all60), condition60,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=all60, colData=coldata, design=~sex+condition60)
rowRanges(ddsSex) <- GR
ddsSex <- DESeq(ddsSex)
resGR10 <- results(ddsSex, format = "GRanges")
resGR10$padj<--log10(resGR10$padj)
d10 <- DataTrack(resGR10, data = "padj", type = "p", name = "-log10(padj)", strand = "*",background.title = "white",col.axis=1,col.title=1,col="black",col.axis=1,col.title=1,ylim=c(0,25))

#Set active scaffolds and then withdraw p values, repeat for next scaffold
chromosome(d21) <- 'Hmel218002o' 
chromosome(d9) <- 'Hmel218002o'
chromosome(d10) <- 'Hmel218002o'
p_raw_d21_Hmel218002o<-values(d21)
p_raw_d9_Hmel218002o<-values(d9)
p_raw_d10_Hmel218002o<-values(d10)

chromosome(d21) <- 'Hmel218003o' 
chromosome(d9) <- 'Hmel218003o' 
chromosome(d10) <- 'Hmel218003o' 
p_raw_d21_Hmel218003o<-values(d21)
p_raw_d9_Hmel218003o<-values(d9)
p_raw_d10_Hmel218003o<-values(d10)

#SET RANGE
write.table(d21@range,"d21.txt")
#Read in ranges and change all factor columns to character columns.
d21_seq<-read.table("d21.txt",header=T)
for(i in 1:length(d21_seq[1,])){
  if(is.factor(d21_seq[,i])){
    d21_seq[,i]<-levels(d21_seq[,i])[d21_seq[,i]]
  }
}
#
write.table(d9@range,"d9.txt")
d9_seq<-read.table("d9.txt",header=T)
for(i in 1:length(d9_seq[1,])){
  if(is.factor(d9_seq[,i])){
    d9_seq[,i]<-levels(d9_seq[,i])[d9_seq[,i]]
  }
}
#
write.table(d10@range,"d10.txt")
d10_seq<-read.table("d10.txt",header=T)
for(i in 1:length(d10_seq[1,])){
  if(is.factor(d10_seq[,i])){
    d10_seq[,i]<-levels(d10_seq[,i])[d10_seq[,i]]
  }
}

#Retrieve chromosome names
seqnames_d21<-d21_seq$seqnames
seqnames_d9<-d9_seq$seqnames
seqnames_d10<-d10_seq$seqnames

#Retrieve mean position of scaffold (mean of beginning and end of scaffold)
mean_pos_d21<-(d21_seq$start+d21_seq$end)/2
mean_pos_d9<-(d9_seq$start+d9_seq$end)/2
mean_pos_d10<-(d10_seq$start+d10_seq$end)/2

#Extract those positions that belong to the two scaffolds you're interested in.
mean_pos_d21_Hmel218002o<-mean_pos_d21[seqnames_d21=="Hmel218002o"]
mean_pos_d9_Hmel218002o<-mean_pos_d9[seqnames_d9=="Hmel218002o"]
mean_pos_d10_Hmel218002o<-mean_pos_d10[seqnames_d10=="Hmel218002o"]

mean_pos_d21_Hmel218003o<-mean_pos_d21[seqnames_d21=="Hmel218003o"]
mean_pos_d9_Hmel218003o<-mean_pos_d9[seqnames_d9=="Hmel218003o"]
mean_pos_d10_Hmel218003o<-mean_pos_d10[seqnames_d10=="Hmel218003o"]

#Fold change
fold_d21_Hmel218002o<--(resGR21$log2FoldChange[seqnames_d=="Hmel218002o"]) #reverse sign of the fold change
fold_d9_Hmel218002o<-resGR9$log2FoldChange[seqnames_d=="Hmel218002o"]
fold_d10_Hmel218002o<-resGR10$log2FoldChange[seqnames_d=="Hmel218002o"]

fold_d21_Hmel218003o<--(resGR21$log2FoldChange[seqnames_d=="Hmel218003o"])
fold_d9_Hmel218003o<-resGR9$log2FoldChange[seqnames_d=="Hmel218003o"]
fold_d10_Hmel218003o<-resGR10$log2FoldChange[seqnames_d=="Hmel218003o"]

#Which positions are free of NA
no_NA_d21_Hmel218002o<-sapply(1:length(fold_d21_Hmel218002o),function(x) sum(is.na(p_raw_d21_Hmel218002o[x]),is.na(mean_pos_d21_Hmel218002o[x]),is.na(fold_d21_Hmel218002o[x])))
no_NA_d9_Hmel218002o<-sapply(1:length(fold_d9_Hmel218002o),function(x) sum(is.na(p_raw_d9_Hmel218002o[x]),is.na(mean_pos_d9_Hmel218002o[x]),is.na(fold_d9_Hmel218002o[x])))
no_NA_d10_Hmel218002o<-sapply(1:length(fold_d10_Hmel218002o),function(x) sum(is.na(p_raw_d10_Hmel218002o[x]),is.na(mean_pos_d10_Hmel218002o[x]),is.na(fold_d10_Hmel218002o[x])))

no_NA_d21_Hmel218003o<-sapply(1:length(fold_d21_Hmel218003o),function(x) sum(is.na(p_raw_d21_Hmel218003o[x]),is.na(mean_pos_d21_Hmel218003o[x]),is.na(fold_d21_Hmel218003o[x])))
no_NA_d9_Hmel218003o<-sapply(1:length(fold_d9_Hmel218003o),function(x) sum(is.na(p_raw_d9_Hmel218003o[x]),is.na(mean_pos_d9_Hmel218003o[x]),is.na(fold_d9_Hmel218003o[x])))
no_NA_d10_Hmel218003o<-sapply(1:length(fold_d10_Hmel218003o),function(x) sum(is.na(p_raw_d10_Hmel218003o[x]),is.na(mean_pos_d10_Hmel218003o[x]),is.na(fold_d10_Hmel218003o[x])))

#Color vectors
col_d21_Hmel218002o_alt<-sapply((1:length(fold_d21_Hmel218002o))[no_NA_d21_Hmel218002o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d21_Hmel218002o[x]<=1.3|(fold_d21_Hmel218002o[x]>=(-1)&fold_d21_Hmel218002o[x]<=1),p_raw_d21_Hmel218002o[x]>1.3&fold_d21_Hmel218002o[x]<(-1),p_raw_d21_Hmel218002o[x]>1.3&fold_d21_Hmel218002o[x]>1)])
col_d9_Hmel218002o_alt<-sapply((1:length(fold_d9_Hmel218002o))[no_NA_d9_Hmel218002o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d9_Hmel218002o[x]<=1.3|(fold_d9_Hmel218002o[x]>=(-1)&fold_d9_Hmel218002o[x]<=1),p_raw_d9_Hmel218002o[x]>1.3&fold_d9_Hmel218002o[x]<(-1),p_raw_d9_Hmel218002o[x]>1.3&fold_d9_Hmel218002o[x]>1)])
col_d10_Hmel218002o_alt<-sapply((1:length(fold_d10_Hmel218002o))[no_NA_d10_Hmel218002o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d10_Hmel218002o[x]<=1.3|(fold_d10_Hmel218002o[x]>=(-1)&fold_d10_Hmel218002o[x]<=1),p_raw_d10_Hmel218002o[x]>1.3&fold_d10_Hmel218002o[x]<(-1),p_raw_d10_Hmel218002o[x]>1.3&fold_d10_Hmel218002o[x]>1)])

col_d21_Hmel218003o_alt<-sapply((1:length(fold_d21_Hmel218003o))[no_NA_d21_Hmel218003o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d21_Hmel218003o[x]<=1.3|(fold_d21_Hmel218003o[x]>=(-1)&fold_d21_Hmel218003o[x]<=1),p_raw_d21_Hmel218003o[x]>1.3&fold_d21_Hmel218003o[x]<(-1),p_raw_d21_Hmel218003o[x]>1.3&fold_d21_Hmel218003o[x]>1)])
col_d9_Hmel218003o_alt<-sapply((1:length(fold_d9_Hmel218003o))[no_NA_d9_Hmel218003o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d9_Hmel218003o[x]<=1.3|(fold_d9_Hmel218003o[x]>=(-1)&fold_d9_Hmel218003o[x]<=1),p_raw_d9_Hmel218003o[x]>1.3&fold_d9_Hmel218003o[x]<(-1),p_raw_d9_Hmel218003o[x]>1.3&fold_d9_Hmel218003o[x]>1)])
col_d10_Hmel218003o_alt<-sapply((1:length(fold_d10_Hmel218003o))[no_NA_d10_Hmel218003o==0],function(x) c("gray89","dodgerblue3","orange")[c(p_raw_d10_Hmel218003o[x]<=1.3|(fold_d10_Hmel218003o[x]>=(-1)&fold_d10_Hmel218003o[x]<=1),p_raw_d10_Hmel218003o[x]>1.3&fold_d10_Hmel218003o[x]<(-1),p_raw_d10_Hmel218003o[x]>1.3&fold_d10_Hmel218003o[x]>1)])


###EFFECT SIZES OVER CHROMOSOME
#Set xlimits
xlims<-c(0,2393341)

#Set ylimits (with maximum y being the rounded up maximum value of all 2*3 plots)
ylims<-c(min(fold_d_Hmel218002o[no_NA_d_Hmel218002o==0],
             fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0],
             fold_d4_Hmel218002o[no_NA_d4_Hmel218002o==0],
             fold_d21_Hmel218002o[no_NA_d21_Hmel218002o==0],
             fold_d9_Hmel218002o[no_NA_d9_Hmel218002o==0],
             fold_d10_Hmel218002o[no_NA_d10_Hmel218002o==0],
             fold_d_Hmel218003o[no_NA_d_Hmel218003o==0],
             fold_d3_Hmel218003o[no_NA_d3_Hmel218003o==0],
             fold_d4_Hmel218003o[no_NA_d4_Hmel218003o==0],
             fold_d21_Hmel218003o[no_NA_d21_Hmel218003o==0],
             fold_d9_Hmel218003o[no_NA_d9_Hmel218003o==0],
             fold_d10_Hmel218003o[no_NA_d10_Hmel218003o==0])
         ,max(fold_d_Hmel218002o[no_NA_d_Hmel218002o==0],
              fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0],
              fold_d4_Hmel218002o[no_NA_d4_Hmel218002o==0],
              fold_d21_Hmel218002o[no_NA_d21_Hmel218002o==0],
              fold_d9_Hmel218002o[no_NA_d9_Hmel218002o==0],
              fold_d10_Hmel218002o[no_NA_d10_Hmel218002o==0],
              fold_d_Hmel218003o[no_NA_d_Hmel218003o==0],
              fold_d3_Hmel218003o[no_NA_d3_Hmel218003o==0],
              fold_d4_Hmel218003o[no_NA_d4_Hmel218003o==0],
              fold_d21_Hmel218003o[no_NA_d21_Hmel218003o==0],
              fold_d9_Hmel218003o[no_NA_d9_Hmel218003o==0],
              fold_d10_Hmel218003o[no_NA_d10_Hmel218003o==0]))


#open pdf where graph will be stored                                
pdf("Hmel218002o_effect_size.pdf")
#Layout to get multiple plots
layout(matrix(1:6,ncol=1))
#spacing around one plot
par(mar=c(2,5,0,0))
#spacing around all three plots
par(oma=c(0,0,1,0.7))
#Create empty plot.
#Give some space on top, on bottom and on left side so that points don't get cut off
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
#Add threshold line
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
#Add points and color if bigger than 1.3 (logarithm in base ten of the p value adjusted)
points(mean_pos_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt=="black"],fold_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218002o_alt[col_d_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt!="black"],fold_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218002o_alt[col_d_Hmel218002o_alt!="black"],0.7))

#Add y axis
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

#Repeat for d3 and d4
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt=="black"],fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218002o_alt[col_d3_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt!="black"],fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218002o_alt[col_d3_Hmel218002o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d4_Hmel218002o[no_NA_d4_Hmel218002o==0][col_d4_Hmel218002o_alt=="black"],fold_d4_Hmel218002o[no_NA_d4_Hmel218002o==0][col_d4_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d4_Hmel218002o_alt[col_d4_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d4_Hmel218002o[no_NA_d4_Hmel218002o==0][col_d4_Hmel218002o_alt!="black"],fold_d4_Hmel218002o[no_NA_d4_Hmel218002o==0][col_d4_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d4_Hmel218002o_alt[col_d4_Hmel218002o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

#Repeat for d21 and d9 and d10
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d21_Hmel218002o[no_NA_d21_Hmel218002o==0][col_d21_Hmel218002o_alt=="black"],fold_d21_Hmel218002o[no_NA_d21_Hmel218002o==0][col_d21_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d21_Hmel218002o_alt[col_d21_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d21_Hmel218002o[no_NA_d21_Hmel218002o==0][col_d21_Hmel218002o_alt!="black"],fold_d21_Hmel218002o[no_NA_d21_Hmel218002o==0][col_d21_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d21_Hmel218002o_alt[col_d21_Hmel218002o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d9_Hmel218002o[no_NA_d9_Hmel218002o==0][col_d9_Hmel218002o_alt=="black"],fold_d9_Hmel218002o[no_NA_d9_Hmel218002o==0][col_d9_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d9_Hmel218002o_alt[col_d9_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d9_Hmel218002o[no_NA_d9_Hmel218002o==0][col_d9_Hmel218002o_alt!="black"],fold_d9_Hmel218002o[no_NA_d9_Hmel218002o==0][col_d9_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d9_Hmel218002o_alt[col_d9_Hmel218002o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d10_Hmel218002o[no_NA_d10_Hmel218002o==0][col_d10_Hmel218002o_alt=="black"],fold_d10_Hmel218002o[no_NA_d10_Hmel218002o==0][col_d10_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d10_Hmel218002o_alt[col_d10_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d10_Hmel218002o[no_NA_d10_Hmel218002o==0][col_d10_Hmel218002o_alt!="black"],fold_d10_Hmel218002o[no_NA_d10_Hmel218002o==0][col_d10_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d10_Hmel218002o_alt[col_d10_Hmel218002o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

#Give only one y axis title
mtext("Fold Change",2,font=2,cex=1.2,line=(-1.9),outer=T)

dev.off()



#scaffold 3
pdf("Hmel218003o_effect_size.pdf")
#### FIRST PLOT (First chromosome)
#Layout to get multiple plots
layout(matrix(1:6,ncol=1))
#spacing around one plot
par(mar=c(2,5,0,0))
#spacing around all three plots
par(oma=c(0,0,1,0.7))
#Create empty plot.
#Give some space on top, on bottom and on left side so that points don't get cut off
#(which they, by the way, do in plotTracks)
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
#Add threshold line
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
#Add points and color if bigger than 1.3
points(mean_pos_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt=="black"],fold_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218003o_alt[col_d_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt!="black"],fold_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218003o_alt[col_d_Hmel218003o_alt!="black"],0.7))
#Add y axis
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

#Repeat for d3 and d4
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt=="black"],fold_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218003o_alt[col_d3_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt!="black"],fold_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218003o_alt[col_d3_Hmel218003o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d4_Hmel218003o[no_NA_d4_Hmel218003o==0][col_d4_Hmel218003o_alt=="black"],fold_d4_Hmel218003o[no_NA_d4_Hmel218003o==0][col_d4_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d4_Hmel218003o_alt[col_d4_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d4_Hmel218003o[no_NA_d4_Hmel218003o==0][col_d4_Hmel218003o_alt!="black"],fold_d4_Hmel218003o[no_NA_d4_Hmel218003o==0][col_d4_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d4_Hmel218003o_alt[col_d4_Hmel218003o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

#Repeat for d21 and d9 and d10
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d21_Hmel218003o[no_NA_d21_Hmel218003o==0][col_d21_Hmel218003o_alt=="black"],fold_d21_Hmel218003o[no_NA_d21_Hmel218003o==0][col_d21_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d21_Hmel218003o_alt[col_d21_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d21_Hmel218003o[no_NA_d21_Hmel218003o==0][col_d21_Hmel218003o_alt!="black"],fold_d21_Hmel218003o[no_NA_d21_Hmel218003o==0][col_d21_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d21_Hmel218003o_alt[col_d21_Hmel218003o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d9_Hmel218003o[no_NA_d9_Hmel218003o==0][col_d9_Hmel218003o_alt=="black"],fold_d9_Hmel218003o[no_NA_d9_Hmel218003o==0][col_d9_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d9_Hmel218003o_alt[col_d9_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d9_Hmel218003o[no_NA_d9_Hmel218003o==0][col_d9_Hmel218003o_alt!="black"],fold_d9_Hmel218003o[no_NA_d9_Hmel218003o==0][col_d9_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d9_Hmel218003o_alt[col_d9_Hmel218003o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),2393341),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(2393341,2393341,2393341),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d10_Hmel218003o[no_NA_d10_Hmel218003o==0][col_d10_Hmel218003o_alt=="black"],fold_d10_Hmel218003o[no_NA_d10_Hmel218003o==0][col_d10_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d10_Hmel218003o_alt[col_d10_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d10_Hmel218003o[no_NA_d10_Hmel218003o==0][col_d10_Hmel218003o_alt!="black"],fold_d10_Hmel218003o[no_NA_d10_Hmel218003o==0][col_d10_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d10_Hmel218003o_alt[col_d10_Hmel218003o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)

#
mtext("Fold Change",2,font=2,cex=1.2,line=(-1.9),outer=T)
dev.off()

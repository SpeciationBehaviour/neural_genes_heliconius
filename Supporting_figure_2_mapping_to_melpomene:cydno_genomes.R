setwd("/Users/matteo/Desktop/")
library("Gviz")
library("DESeq2")
options(ucscChromosomeNames=FALSE)


#make genomic ranges
annotation<-read.table("gene_info.txt",header=TRUE) 
annotation <- annotation[order(annotation$gene_id),]
GR<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))
a <- AnnotationTrack(GR, name = "gene ranges") 

#
annotation_cyd<-read.table("gene_info_cydno.txt",header=TRUE)
annotation_cyd <- annotation_cyd[order(annotation_cyd$gene_id),]
GR2<-makeGRangesFromDataFrame(annotation_cyd, seqnames.field=c("scaffold"))
a3 <- AnnotationTrack(GR2, name = "gene ranges") 




#################################DIFFERENTIAL EXPRESSION

########################Adults TO MELPOMENE
allAdult<-read.table("Adults_toMP_gene_counts.txt",header=TRUE, row.names = 1)
allAdult<-allAdult[,c(1:23)]  
allAdult<- as.matrix(allAdult)
(sex <- factor(c(rep("male",7), rep("female", 5),rep("male",4), rep("female",7))))
(conditionAdult <- factor(c(rep("mel", 12), rep("cydno", 11))))
(coldata <- data.frame(row.names=colnames(allAdult), conditionAdult,sex))
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+conditionAdult)
rowRanges(ddsSex) <- GR #add genomic ranges to differential expression results
dds <- DESeq(ddsSex)
resGR <- results(dds, format = "GRanges")
resGR$padj<--log10(resGR$padj) #transform p values adjusted in -log base10
d <- DataTrack(resGR, data = "padj", type = "p", name = "-log10(padj)", strand = "*") #store p-values and genomic ranges in Gvizobject


########################Adults TO CYDNO
allAdult2<-read.table("Adults_toCP_gene_counts.txt",header=TRUE, row.names = 1)
allAdult2<-allAdult2[,c(1:23)]  
allAdult2<- as.matrix(allAdult2)
(sex <- factor(c(rep("male",7), rep("female", 5),rep("male",4), rep("female",7))))
(conditionAdult2 <- factor(c(rep("mel", 12), rep("cydno", 11))))
(coldata2 <- data.frame(row.names=colnames(allAdult2), conditionAdult2,sex))
ddsSex2 <- DESeqDataSetFromMatrix(countData=allAdult2, colData=coldata2, design=~sex+conditionAdult2)
rowRanges(ddsSex2) <- GR2 #add genomic ranges to differential expression results
dds2 <- DESeq(ddsSex2)
resGR3 <- results(dds2, format = "GRanges")
resGR3$padj<--log10(resGR3$padj) #transform p values adjusted in -log base10
d3 <- DataTrack(resGR3, data = "padj", type = "p", name = "-log10(padj)", strand = "*") #store p-values and genomic ranges in Gvizobject



#######################

#Set active Chromosome and then withdraw p values
chromosome(d) <- 'Hmel218002o' 
chromosome(d3) <- 'Hmel218002o' 
p_raw_d_Hmel218002o<-values(d)
p_raw_d3_Hmel218002o<-values(d3)

chromosome(d) <- 'Hmel218003o' 
chromosome(d3) <- 'Hmel218003o' 
p_raw_d_Hmel218003o<-values(d)
p_raw_d3_Hmel218003o<-values(d3)

#Safe ranged
write.table(d@range,"d.txt")
#Red in ranges and change all factor columns to character columns.
d_seq<-read.table("d.txt",header=T)
for(i in 1:length(d_seq[1,])){
  if(is.factor(d_seq[,i])){
    d_seq[,i]<-levels(d_seq[,i])[d_seq[,i]]
  }
}

#Repeat the same for other two
write.table(d3@range,"d3.txt")
d3_seq<-read.table("d3.txt",header=T)
for(i in 1:length(d3_seq[1,])){
  if(is.factor(d3_seq[,i])){
    d3_seq[,i]<-levels(d3_seq[,i])[d3_seq[,i]]
  }
}     

#Retrieve chromosome names
seqnames_d<-d_seq$seqnames
seqnames_d3<-d3_seq$seqnames

#Retrieve mean position of scaffold (mean of beginning and end of scaffold)
mean_pos_d<-(d_seq$start+d_seq$end)/2
mean_pos_d3<-(d3_seq$start+d3_seq$end)/2      


#Extract those positions that belong to the two chromosomes you're interested in.
mean_pos_d_Hmel218002o<-mean_pos_d[seqnames_d=="Hmel218002o"]
mean_pos_d3_Hmel218002o<-mean_pos_d3[seqnames_d3=="Hmel218002o"]

mean_pos_d_Hmel218003o<-mean_pos_d[seqnames_d=="Hmel218003o"]
mean_pos_d3_Hmel218003o<-mean_pos_d3[seqnames_d3=="Hmel218003o"]


#Fold change
fold_d_Hmel218002o<-resGR$log2FoldChange[seqnames_d=="Hmel218002o"]
fold_d3_Hmel218002o<-resGR3$log2FoldChange[seqnames_d3=="Hmel218002o"]

fold_d_Hmel218003o<-resGR$log2FoldChange[seqnames_d=="Hmel218003o"]
fold_d3_Hmel218003o<-resGR3$log2FoldChange[seqnames_d3=="Hmel218003o"]


#Which positions are free of NA
no_NA_d_Hmel218002o<-sapply(1:length(fold_d_Hmel218002o),function(x) sum(is.na(p_raw_d_Hmel218002o[x]),is.na(mean_pos_d_Hmel218002o[x]),is.na(fold_d_Hmel218002o[x])))
no_NA_d3_Hmel218002o<-sapply(1:length(fold_d3_Hmel218002o),function(x) sum(is.na(p_raw_d3_Hmel218002o[x]),is.na(mean_pos_d3_Hmel218002o[x]),is.na(fold_d3_Hmel218002o[x])))

no_NA_d_Hmel218003o<-sapply(1:length(fold_d_Hmel218003o),function(x) sum(is.na(p_raw_d_Hmel218003o[x]),is.na(mean_pos_d_Hmel218003o[x]),is.na(fold_d_Hmel218003o[x])))
no_NA_d3_Hmel218003o<-sapply(1:length(fold_d3_Hmel218003o),function(x) sum(is.na(p_raw_d3_Hmel218003o[x]),is.na(mean_pos_d3_Hmel218003o[x]),is.na(fold_d3_Hmel218003o[x])))


#Color vectors
col_d_Hmel218002o_alt<-sapply((1:length(fold_d_Hmel218002o))[no_NA_d_Hmel218002o==0],function(x) c("black","dodgerblue","gold")[c(p_raw_d_Hmel218002o[x]<=1.3|(fold_d_Hmel218002o[x]>=(-1)&fold_d_Hmel218002o[x]<=1),p_raw_d_Hmel218002o[x]>1.3&fold_d_Hmel218002o[x]<(-1),p_raw_d_Hmel218002o[x]>1.3&fold_d_Hmel218002o[x]>1)])
col_d3_Hmel218002o_alt<-sapply((1:length(fold_d3_Hmel218002o))[no_NA_d3_Hmel218002o==0],function(x) c("black","dodgerblue","gold")[c(p_raw_d3_Hmel218002o[x]<=1.3|(fold_d3_Hmel218002o[x]>=(-1)&fold_d3_Hmel218002o[x]<=1),p_raw_d3_Hmel218002o[x]>1.3&fold_d3_Hmel218002o[x]<(-1),p_raw_d3_Hmel218002o[x]>1.3&fold_d3_Hmel218002o[x]>1)])

col_d_Hmel218003o_alt<-sapply((1:length(fold_d_Hmel218003o))[no_NA_d_Hmel218003o==0],function(x) c("black","dodgerblue","gold")[c(p_raw_d_Hmel218003o[x]<=1.3|(fold_d_Hmel218003o[x]>=(-1)&fold_d_Hmel218003o[x]<=1),p_raw_d_Hmel218003o[x]>1.3&fold_d_Hmel218003o[x]<(-1),p_raw_d_Hmel218003o[x]>1.3&fold_d_Hmel218003o[x]>1)])
col_d3_Hmel218003o_alt<-sapply((1:length(fold_d3_Hmel218003o))[no_NA_d3_Hmel218003o==0],function(x) c("black","dodgerblue","gold")[c(p_raw_d3_Hmel218003o[x]<=1.3|(fold_d3_Hmel218003o[x]>=(-1)&fold_d3_Hmel218003o[x]<=1),p_raw_d3_Hmel218003o[x]>1.3&fold_d3_Hmel218003o[x]<(-1),p_raw_d3_Hmel218003o[x]>1.3&fold_d3_Hmel218003o[x]>1)])


#Set xlimits
xlims<-c(0,1200000)

#Set xlimits
#xlims<-c(0,2393341)

#Set ylimits (with maximum y being the rounded up maximum value of all 2*3 plots)
#If you really want the y area to end at 30, just set this to c(0,30)
ylims<-c(min(fold_d_Hmel218002o[no_NA_d_Hmel218002o==0],
             fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0])
         ,max(fold_d_Hmel218002o[no_NA_d_Hmel218002o==0],
              fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0]))




############## SCAFFOLD 2
pdf("Hmel218002o_ALL.pdf")
plotTracks(list(d,d3), chromosome = "Hmel218002o", from = 0,to = 1200000, ylim=c(0,35))

#### FIRST PLOT (First chromosome)
#Layut to get multiple plots
layout(matrix(1:2,ncol=1))
#spacing around one plot
par(mar=c(2,5,0,0))
#spacing around all three plots
par(oma=c(0,0,1,0.7))


#ADULTS #to MELPOMENE
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),1200000),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
#Add threshold line
segments(c(0,0,0),c(-1,1,0),c(1200000,1200000,1200000),c(-1,1,0),lty=c("dotted","dotted","solid"))
#Add points and color if bigger than 1.3
points(mean_pos_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt=="black"],fold_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218002o_alt[col_d_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt!="black"],fold_d_Hmel218002o[no_NA_d_Hmel218002o==0][col_d_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218002o_alt[col_d_Hmel218002o_alt!="black"],0.7))

#Add y axis
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)


#ADULTS #to CYDNO
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),1200000),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(1200000,1200000,1200000),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt=="black"],fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218002o_alt[col_d3_Hmel218002o_alt=="black"],0.7))
points(mean_pos_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt!="black"],fold_d3_Hmel218002o[no_NA_d3_Hmel218002o==0][col_d3_Hmel218002o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218002o_alt[col_d3_Hmel218002o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)


dev.off()





#########################GRAPH SCAFFOLD 3

pdf("Hmel218003o_ALL.pdf")
plotTracks(list(d,d3), chromosome = "Hmel218003o", from = 0,to = 1200000, ylim=c(0,35))

#### FIRST PLOT (First chromosome)
#Layut to get multiple plots
layout(matrix(1:2,ncol=1))
#spacing around one plot
par(mar=c(2,5,0,0))
#spacing around all three plots
par(oma=c(0,0,1,0.7))


#MELPOMENE
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),1200000),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
#Add threshold line
segments(c(0,0,0),c(-1,1,0),c(1200000,1200000,1200000),c(-1,1,0),lty=c("dotted","dotted","solid"))
#Add points and color if bigger than 1.3
points(mean_pos_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt=="black"],fold_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218003o_alt[col_d_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt!="black"],fold_d_Hmel218003o[no_NA_d_Hmel218003o==0][col_d_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d_Hmel218003o_alt[col_d_Hmel218003o_alt!="black"],0.7))
#Add y axis
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)


##CYDNO
plot(1,1,
     xlim=c(xlims[1]-0.01*(abs(xlims[1]-xlims[2])),1200000),
     ylim=c(ylims[1]-0.02*(abs(ylims[1]-ylims[2])),
            ylims[2]+0.01*(abs(ylims[1]-ylims[2]))),
     type="n",xlab="",ylab="",xaxt="none",yaxt="none",axes=F,xaxs="i",yaxs="i")
segments(c(0,0,0),c(-1,1,0),c(1200000,1200000,1200000),c(-1,1,0),lty=c("dotted","dotted","solid"))
points(mean_pos_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt=="black"],fold_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt=="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218003o_alt[col_d3_Hmel218003o_alt=="black"],0.7))
points(mean_pos_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt!="black"],fold_d3_Hmel218003o[no_NA_d3_Hmel218003o==0][col_d3_Hmel218003o_alt!="black"],
       pch=21,lwd=0.5,cex=1.5,bg=adjustcolor(col_d3_Hmel218003o_alt[col_d3_Hmel218003o_alt!="black"],0.7))
axis(2,line=0.2,las=1,lwd=1.5,font=2)
axis(2,line=0.2,lwd=1.5,at=ylims,tck=0,labels=F)


dev.off()






######see which GENES diff epressed in both comparisons
ddsSex <- DESeqDataSetFromMatrix(countData=allAdult, colData=coldata, design=~sex+conditionAdult)
dds <- DESeq(ddsSex)
res <- results(dds)
resSig <- subset(res, padj < 0.05) 
#write.table(resSig, file="Adults_DE_toMP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diffAdult_to_mel<-read.table("Adults_DE_toMP.txt",header=TRUE)
upmel<-diffAdult_to_mel[diffAdult_to_mel$log2FoldChange >1,]
upcydno<-diffAdult_to_mel[diffAdult_to_mel$log2FoldChange< (-1),]
all<-rbind(upmel,upcydno)
all <- merge(all,annotation,by="gene_id")   

#
ddsSex2 <- DESeqDataSetFromMatrix(countData=allAdult2, colData=coldata2, design=~sex+conditionAdult2)
dds2 <- DESeq(ddsSex2)
res2 <- results(dds2)
resSig2 <- subset(res2, padj < 0.05) 
#write.table(resSig2, file="Adults_DE_toCP.txt", row.names = TRUE, quote=FALSE)  #add gene_id file first column

diffAdult_to_cyd<-read.table("Adults_DE_toCP.txt",header=TRUE)
upmel2<-diffAdult_to_cyd[diffAdult_to_cyd$log2FoldChange >1,]
upcydno2<-diffAdult_to_cyd[diffAdult_to_cyd$log2FoldChange< (-1),]
all2<-rbind(upmel2,upcydno2)
all2<- merge(all2,annotation_cyd,by="gene_id") 


#genes diff expressed until optix #CHR18
#to MP
only18_1<-all[all$scaffold == "Hmel218001o",] #
only18_2<-all[all$scaffold == "Hmel218002o",] #
only18_3<-all[all$scaffold =="Hmel218003o",] #   
untiloptix<-only18_3[only18_3$end <800000,] #

#to CP
only18_1_toCP<-all2[all2$scaffold == "Hmel218001o",] #
only18_2_toCP<-all2[all2$scaffold == "Hmel218002o",] #
only18_3_toCP<-all2[all2$scaffold =="Hmel218003o",] #   
untiloptix_toCP<-only18_3_toCP[only18_3_toCP$end <800000,] # 


#genes diff expressed in QTL candidate regions ("QTL 1.5 lod score")
#to MP
only18_3<-all[all$scaffold =="Hmel218003o",] #   
untilqtl<-only18_3[only18_3$end <2393341,] # 

#to CP
only18_3_toCP<-all2[all2$scaffold =="Hmel218003o",] #   
untilqtl_toCP<-only18_3_toCP[only18_3_toCP$end <2393341,] #   


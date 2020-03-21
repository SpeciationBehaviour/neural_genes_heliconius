setwd("/Users/matteo/Desktop/")

#
library("Gviz")
options(ucscChromosomeNames=FALSE)

#make genomic ranges MELPOMENE
annotation<-read.table("gene_info.txt",header=TRUE) 
#in dds2 results ordered by gene name, so order gene ranges also by gene name
annotation <- annotation[order(annotation$gene_id),]
GR<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))
a <- AnnotationTrack(GR, name = "gene ranges")


##### SIMON MARTIN'S DATA

########100kb windows
admixture<-read.csv("bar92.DP8MP4BIMAC2HET75.fourPopPol_melG_melW_cyd_num.w100m1s20.merged.csv",header=TRUE) #13557 observations (with NAs)
#where D is negative fd statistics are meaningless = no evidence for admixture  #set to 0! (not to NA)
Dnegative<-admixture[as.character(admixture$D)<0,] 
admixture$D[admixture$D < 0] <- 0
#add genomic range to dds object
GR2<-makeGRangesFromDataFrame(admixture, seqnames.field=c("scaffold"))
Fd<-as.matrix(admixture[,c(12)])
Fd<-data.frame(Fd)
Fd$Fd[Fd$Fd < 0] <- 0   #trasform negative Fd values to zero
values(GR2) <- DataFrame(Fd$Fd)
df100 <- DataTrack(GR2, data = "Fd.Fd", strand = "*",type = c("a"),lwd=2,col="black",background.title = "white",col.axis=1,ylim=c(0,0.9),cex.axis=1)


###############20 kb windows
admixture2<-read.csv("bar92.DP8MP4BIMAC2HET75.fourPopPol_melG_melW_cyd_num.w20m300s5.csv",header=TRUE) 
#add genomic ranges
GR2<-makeGRangesFromDataFrame(admixture2, seqnames.field=c("scaffold"))
Fd<-as.matrix(admixture2[,c(12)])
Fd<-data.frame(Fd)
Fd$Fd[Fd$Fd < 0] <- 0   #trasform negative Fd values to zero
values(GR2) <- DataFrame(Fd$Fd)
df20 <- DataTrack(GR2, data = "Fd.Fd", strand = "*",type = c("a"),lwd=2,col="black",background.title = "white",col.axis=1,ylim=c(0,0.9),cex.axis=1)



########################GRAPHs
 
#mel- cydno - 20kb
plotTracks(list(df20), chromosome = "chr18", from = 0,to = 2746316) 
#WHOLE CHROMOSOME #with 100kb
plotTracks(list(df100), chromosome = "chr18") 
#WHOLE CHROMOSOME #with both
plotTracks(list(df100, df20), chromosome = "chr18", from = 0,to = 2746316) 

#to see where boundaries
plotTracks(list(df100), chromosome = "chr18", from = 0,to = 628045) #peak
plotTracks(list(df100), chromosome = "chr18", from = 0,to = 2746316) #interval
plotTracks(list(df20), chromosome = "chr18", from = 0,to = 628045) #peak
plotTracks(list(df20), chromosome = "chr18", from = 0,to = 1000000) #1Mb
plotTracks(list(df20), chromosome = "chr18", from = 0,to = 2000000) #2Mb
#to see where genes
annotation<-read.table("gene_info.txt",header=TRUE) 

######
gene2<-annotation
gene2<-gene2[gene2$scaffold == "Hmel218002o",]
gene2$scaffold<- "chr18"
gene2$start<-gene2$start + 47928
gene2$end<-gene2$end + 47928

gene3<-annotation
gene3<-gene3[gene3$scaffold == "Hmel218003o",]
gene3$scaffold<- "chr18"
gene3$start<-gene3$start + 352975
gene3$end<-gene3$end + 352975

annotation<-rbind(gene2,gene3)
annotation <- annotation[order(annotation$gene_id),]
GR<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))
a <- AnnotationTrack(GR, name = "gene ranges")
#can plot with "a" track #annotation

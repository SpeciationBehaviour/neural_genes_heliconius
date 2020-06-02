setwd("/Users/matteo/Desktop/")

library("Gviz")
options(ucscChromosomeNames=FALSE)

#make genomic ranges to match gene physical position and estimated admixture proportions
annotation<-read.table("gene_info.txt",header=TRUE) 
#in dds2 results ordered by gene name, so order gene ranges also by gene name
annotation <- annotation[order(annotation$gene_id),]
GR<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))
a <- AnnotationTrack(GR, name = "gene ranges")

### SIMON MARTIN'S ADMIXTURE PROPORTIONS
#Estimated in 100kb windows:
admixture<-read.csv("bar92.DP8MP4BIMAC2HET75.fourPopPol_melG_melW_cyd_num.w100m1s20.merged.csv",header=TRUE) #13557 observations (with NAs)
#where D is negative fd statistics are meaningless = no evidence for admixture  #set to 0! 
Dnegative<-admixture[as.character(admixture$D)<0,] 
admixture$D[admixture$D < 0] <- 0
#add genomic range to dds object
GR2<-makeGRangesFromDataFrame(admixture, seqnames.field=c("scaffold"))
Fd<-as.matrix(admixture[,c(12)])
Fd<-data.frame(Fd)
Fd$Fd[Fd$Fd < 0] <- 0   #trasform negative Fd values to zero
values(GR2) <- DataFrame(Fd$Fd)
#store values in DataTRack object
df100 <- DataTrack(GR2, data = "Fd.Fd", strand = "*",type = c("a"),lwd=2,col="black",background.title = "white",col.axis=1,ylim=c(0,0.9),cex.axis=1)

#Estimated in 20kb windows:
admixture2<-read.csv("bar92.DP8MP4BIMAC2HET75.fourPopPol_melG_melW_cyd_num.w20m300s5.csv",header=TRUE) 
#add genomic ranges
GR2<-makeGRangesFromDataFrame(admixture2, seqnames.field=c("scaffold"))
Fd<-as.matrix(admixture2[,c(12)])
Fd<-data.frame(Fd)
Fd$Fd[Fd$Fd < 0] <- 0   #trasform negative Fd values to zero
values(GR2) <- DataFrame(Fd$Fd)
df20 <- DataTrack(GR2, data = "Fd.Fd", strand = "*",type = c("a"),lwd=2,col="black",background.title = "white",col.axis=1,ylim=c(0,0.9),cex.axis=1)

#PRODUCE GRAPHs
#QTL #with admixture proportions in 20kb windows
plotTracks(list(df20), chromosome = "chr18", from = 0,to = 2746316) 
#WHOLE CHROMOSOME #with admixture proportions in 100kb windows
plotTracks(list(df100), chromosome = "chr18") 

#In case one wants to overlay gene models below the graph
annotation<-read.table("gene_info.txt",header=TRUE) 

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
#can plot with "a" track (graph is informative only for small ranges..)


#to make fd distribution graph
admixture2<-read.csv("bar92.DP8MP4BIMAC2HET75.fourPopPol_melG_melW_cyd_num.w20m300s5.csv",header=TRUE) 
admixture2 <- admixture2[admixture2$fd >= 0, ]
distribution<-density(admixture2$fd, na.rm=TRUE)
plot(distribution, lwd=2)
#fd values
#fd values > outmost fd value> within 20 kb upstream #+center of protein-coding part
abline (v=0.0380,lty=2,col="black") #at 1038778
abline (v=0.0032,lty=2,col="black") # optix  #HMEL001028 chr18 1058778 1059581
abline (v=0.0017,lty=2,col="black") #at 1062581

abline (v=0.0108,lty=2,col="goldenrod1") #314856
abline (v=0.0126,lty=2,col="goldenrod1") #HMEL009992g4 chr18 334856 345144
abline (v=0.1259,lty=2,col="goldenrod1") #365144

abline (v=0.1955,lty=2,col="goldenrod1") # 453901
abline (v=0.1331,lty=2,col="goldenrod1") #regucalcin2 #HMEL013552g1 chr18 473901 475137 +  HMEL034199g1 chr18 478430 479866
abline (v=0.0515,lty=2,col="goldenrod1") #495137

abline (v=0.2250,lty=2,col="dodgerblue3") #448445
abline (v=0.1405,lty=2,col="dodgerblue3") #regucalcin1 #HMEL013551g4 chr18 468445 472652
abline (v=0.0515,lty=2,col="dodgerblue3") #492652

abline (v=0.1655,lty=2,col="dodgerblue3") #435434
abline (v=0.2250,lty=2,col="dodgerblue3") #WD40 #HMEL013551g3 chr18 455434 459952
abline (v=0.1171,lty=2,col="dodgerblue3") #479952

abline (v=0.0185,lty=2,col="dodgerblue3") #136919
abline (v=0.1530,lty=2,col="dodgerblue3") #cyst protease #HMEL009684g2 chr18 156919 162617
abline (v=0.0267,lty=2,col="dodgerblue3")  #182617

abline (v=0.0620,lty=2,col="dodgerblue3") # at 2700001 2720000
abline (v=0.3184,lty=2,col="dodgerblue3") #MORN #chr18 2720420 2730566
abline (v=0.1282,lty=2,col="dodgerblue3")#


#quantiles
admixture2<-admixture2[!is.na(admixture2$fd), ]
quantile(admixture2$fd, probs = c(0.05)) #0.0306 
quantile(admixture2$fd, probs = c(0.25)) #0.1488
abline (v=0.0306,lty=2, col="black")
abline (v=0.1488,lty=2, col="black")

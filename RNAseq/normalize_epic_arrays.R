#### This analysis is for analysis of DNA EPIC array analysis. NOTE!!! There is not an annotation available for hg38, so the genomic coordinates are hg19.

setwd("Z:/Wendy_Kellner/DNMT1/DNMT1_AML_Epic_methylation_NYU")

library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(ggplot2)

######################## Reading in files ###########################
###Read in the sample sheet with the phenotypic data
targets <- read.metharray.sheet(base = "path_to_idat")
RGset <- read.metharray.exp(targets = targets, force = TRUE)
targets$ID <- paste(targets$CellLine,targets$Dose,targets$Treatment,targets$Time,targets$Details,sep=".")
sampleNames(RGset) <- targets$ID

################Annotate data########################################
annotation(RGset)
RGset@annotation=c(array='IlluminaHumanMethylationEPIC', annotation='ilm10b2.hg19')

######################### Remove poor qualtiy samples #################################
detP <- detectionP(RGset)
keep <- colMeans(detP) < 0.05
RGset <- RGset[,keep]

########### Processing and normalization of data ###########################
mSetSq <- preprocessSWAN(RGset)      
MSet.raw <- preprocessFunnorm(RGset) 

########### make the violin plot for publication ###########################
 var<-data.frame(getBeta(mSetSq))
 subs<-var[(1:100000),]
 subs2<-(rep(c("DMSO","GSK762"),7))
 subs<-rbind(subs2,subs)
 row.names(subs[1,])<-c("Treatment")
 df.m <- reshape2::melt(subs, id=subs[1,])
 p<-ggplot(df.m, aes(x = variable, y = value),fill=subs) + geom_violin()
    +scale_fill_manual(values=c("blue","red"))
 p + stat_summary(aes(group=1),fun.y=mean, geom="point", size=10,shape=95,col=c("blue","red","blue","red","blue","red","blue","red","blue","red","blue","red","blue","red"))
       

       
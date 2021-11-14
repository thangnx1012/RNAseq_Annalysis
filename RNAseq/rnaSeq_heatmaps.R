setwd("Z:/Wendy_Kellner/DNMT1/DNMT1_AML/DNMT1_AML_RNAseq_NYU/DNMT1_AML_heatmaps")

require ("gplots")
require ("RColorBrewer")

library(gplots)
library(RColorBrewer)

# import the symbols of your genelist
genes <- read.delim("DNMT1_MV411_d4_GSK_treatments_sig_genes.txt")
####Get the column of the genelist you are interested in making a heatmap###
gene_list <- data.frame(genes$MV411_400nM032_d4)

#import the expression output from limma
gene_exp <- read.delim("DNMT1_MV411_lfc_table.txt")

# get only expression values from your gene list
values <- merge(gene_exp, gene_list, by.x="X",by.y="genes.MV411_400nM032_d4", sort=FALSE)
values <- values[!duplicated(values$X),]
row.names(values)<-values[,1]
values<-values[,-1]

write.table(values,file="MV411_d4_GSK_treatments_lfc_heatmap.txt",sep="\t",quote=FALSE)

###only interested in the 400nM of GSK032 and 1000nM GSK862 from day4. Remove all other columns###
sub <- values[ ,c(7,9,10)]
mat_data <-data.matrix(sub)


my_palette <- colorRampPalette(c("navyblue", "white", "red"))(n = 179)
col_breaks = c(seq(-4.0,-1.0,length=60),  # for blue
              seq(-0.99,0.99,length=60),              # for white
               seq(1.0,4.0,length=60))              # for red

heatmap.2(mat_data,
          main = "DNMT1 AML Fold Change", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=NULL,     # only draw a row dendrogram
          Rowv =TRUE,
          Colv= FALSE,
          key = TRUE,
          lmat = rbind(c(0,3),c(0,4),c(2,1)),
          lwid = c(1.5,4),
          lhei = c(0.5,0.8,8),
          keysize = 0.8,
          cexCol = 2,
          srtCol = 90,
          offsetCol = 0.2)

#####################################################

####Get the column of the genelist you are interested in making a heatmap###
gene_list <- data.frame(genes$MV411_400nM032_DAC_overlap_increased_d4)

# get only expression values from your gene list
values <- merge(gene_exp, gene_list, by.x="X",by.y="genes.MV411_400nM032_DAC_overlap_increased_d4", sort=FALSE)
values <- values[!duplicated(values$X),]
row.names(values)<-values[,1]
values<-values[,-1]

###only interested in the 400nM of GSK032 and 1000nM GSK862 from day4. Remove all other columns###
sub <- values[ ,c(13:24)]
mat_data <-data.matrix(sub)


my_palette <- colorRampPalette(c("navyblue", "white", "red"))(n = 179)
col_breaks = c(seq(-4.0,-1.0,length=60),  # for blue
               seq(-0.99,0.99,length=60),              # for white
               seq(1.0,4.0,length=60))              # for red

heatmap.2(mat_data,
          main = "DNMT1 AML Fold Change", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=NULL,     # only draw a row dendrogram
          Rowv =TRUE,
          Colv= FALSE,
          key = TRUE,
          lmat = rbind(c(0,3),c(0,4),c(2,1)),
          lwid = c(1.5,4),
          lhei = c(0.5,0.8,8),
          keysize = 0.8,
          cexCol = 2,
          srtCol = 90,
          offsetCol = 0.2)

#####################################################

####Get the column of the genelist you are interested in making a heatmap###
gene_list <- data.frame(genes$MV411_400nM032_DAC_overlap_decreased_d4)

# get only expression values from your gene list
values <- merge(gene_exp, gene_list, by.x="X",by.y="genes.MV411_400nM032_DAC_overlap_decreased_d4", sort=FALSE)
values <- values[!duplicated(values$X),]
row.names(values)<-values[,1]
values<-values[,-1]

###only interested in the 400nM of GSK032 and 1000nM GSK862 from day4. Remove all other columns###
sub <- values[ ,c(13:24)]
mat_data <-data.matrix(sub)


my_palette <- colorRampPalette(c("navyblue", "white", "red"))(n = 179)
col_breaks = c(seq(-4.0,-1.0,length=60),  # for blue
               seq(-0.99,0.99,length=60),              # for white
               seq(1.0,4.0,length=60))              # for red

heatmap.2(mat_data,
          main = "DNMT1 AML Fold Change", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=NULL,     # only draw a row dendrogram
          Rowv =TRUE,
          Colv= FALSE,
          key = TRUE,
          lmat = rbind(c(0,3),c(0,4),c(2,1)),
          lwid = c(1.5,4),
          lhei = c(0.5,0.8,8),
          keysize = 0.8,
          cexCol = 2,
          srtCol = 90,
          offsetCol = 0.2)

#####################################################

####Get the column of the genelist you are interested in making a heatmap###
gene_list <- data.frame(genes$MV411_10000nM032_DAC_overlap_increased_d4)

# get only expression values from your gene list
values <- merge(gene_exp, gene_list, by.x="X",by.y="genes.MV411_10000nM032_DAC_overlap_increased_d4", sort=FALSE)
values <- values[!duplicated(values$X),]
row.names(values)<-values[,1]
values<-values[,-1]

###only interested in the 400nM of GSK032 and 1000nM GSK862 from day4. Remove all other columns###
sub <- values[ ,c(13:24)]
mat_data <-data.matrix(sub)


my_palette <- colorRampPalette(c("navyblue", "white", "red"))(n = 179)
col_breaks = c(seq(-4.0,-1.0,length=60),  # for blue
               seq(-0.99,0.99,length=60),              # for white
               seq(1.0,4.0,length=60))              # for red

heatmap.2(mat_data,
          main = "DNMT1 AML Fold Change", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=NULL,     # only draw a row dendrogram
          Rowv =TRUE,
          Colv= FALSE,
          key = TRUE,
          lmat = rbind(c(0,3),c(0,4),c(2,1)),
          lwid = c(1.5,4),
          lhei = c(0.5,0.8,8),
          keysize = 0.8,
          cexCol = 2,
          srtCol = 90,
          offsetCol = 0.2)

#####################################################

####Get the column of the genelist you are interested in making a heatmap###
gene_list <- data.frame(genes$MV411_10000nM032_DAC_overlap_decreased_d4)

# get only expression values from your gene list
values <- merge(gene_exp, gene_list, by.x="X",by.y="genes.MV411_10000nM032_DAC_overlap_decreased_d4", sort=FALSE)
values <- values[!duplicated(values$X),]
row.names(values)<-values[,1]
values<-values[,-1]

###only interested in the 400nM of GSK032 and 1000nM GSK862 from day4. Remove all other columns###
sub <- values[ ,c(13:24)]
mat_data <-data.matrix(sub)


my_palette <- colorRampPalette(c("navyblue", "white", "red"))(n = 179)
col_breaks = c(seq(-4.0,-1.0,length=60),  # for blue
               seq(-0.99,0.99,length=60),              # for white
               seq(1.0,4.0,length=60))              # for red

png(filename="Dose_response" )
heatmap.2(mat_data,
          main = "DNMT1 AML Fold Change", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=NULL,     # only draw a row dendrogram
          Rowv =TRUE,
          Colv= FALSE,
          key = TRUE,
          lmat = rbind(c(0,3),c(0,4),c(2,1)),
          lwid = c(1.5,4),
          lhei = c(0.5,0.8,8),
          keysize = 0.8,
          cexCol = 2,
          srtCol = 90,
          offsetCol = 0.2)
dev.off()

##############################################
####Get the column of the genelist you are interested in making a heatmap###
gene_list <- data.frame(genes$MV411_GSK032_DAC_dose_all_genes)

# get only expression values from your gene list
values <- merge(gene_exp, gene_list, by.x="X",by.y="genes.MV411_GSK032_DAC_dose_all_genes", sort=FALSE)
values <- values[!duplicated(values$X),]
row.names(values)<-values[,1]
values<-values[,-1]

###only interested in the 400nM of GSK032 and 1000nM GSK862 from day4. Remove all other columns###
sub <- values[ ,c(7,9,10)]
mat_data <-data.matrix(sub)


my_palette <- colorRampPalette(c("navyblue", "white", "red"))(n = 179)
col_breaks = c(seq(-4.0,-1.0,length=60),  # for blue
               seq(-0.99,0.99,length=60),              # for white
               seq(1.0,4.0,length=60))              # for red

heatmap.2(mat_data,
          main = "DNMT1 AML Fold Change", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=NULL,     # only draw a row dendrogram
          Rowv =TRUE,
          Colv= FALSE,
          key = TRUE,
          lmat = rbind(c(0,3),c(0,4),c(2,1)),
          lwid = c(1.5,4),
          lhei = c(0.5,0.8,8),
          keysize = 0.8,
          cexCol = 2,
          srtCol = 90,
          offsetCol = 0.2)

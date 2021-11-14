#Example run of DESeq2 to quantify differential expression of ERVs
#Input is a count table (cntTbale) with hERV names, and read counts for 2 control and 2 treated samples
#This code was run for each comparison

data <- read.table(count_table,header=T,row.names=1)
groups <- factor(c(rep("TGroup",2),rep("CGroup",2)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
library(DESeq2, quietly=T)
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$condition = relevel(dds$groups,"CGroup")
dds <- DESeq(dds)
res <- results(dds)
write.table(res, file=outFil1, sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) & (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file=outFile2,sep="\t", quote=F)


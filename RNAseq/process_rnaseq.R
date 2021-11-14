#The input to this script is the output of Salmon. The output will be one directory per sample including a quan.sf file
#Example comparisons o THP1 are included in this script, it can be run for each cell line of interest.
#To generate a command to run salmon for each sample, run the below for loop in bash with the "fastq" directory being where the raw files from this publication are located
#for f in fastq/*L001_R1_001.fastq.gz; do echo "salmon quant -i ~/Desktop/Annotations/Salmon_gencodeGRCh38v23_index/ -l A -1 $f $(echo $f | sed s/_L001_R1/_L002_R1/g) $(echo $f | sed s/_L001_R1/_L003_R1/g) $(echo $f | sed s/_L001_R1/_L004_R1/g) $(echo $f | sed s/_L001_R1/_L005_R1/g) $(echo $f | sed s/_L001_R1/_L006_R1/g) -2 $(echo $f | sed s/L001_R1/L001_R2/g) $(echo $f | sed s/_L001_R1/_L002_R2/g) $(echo $f | sed s/_L001_R1/_L003_R2/g) $(echo $f | sed s/_L001_R1/_L004_R2/g) $(echo $f | sed s/_L001_R1/_L005_R2/g) $(echo $f | sed s/_L001_R1/_L006_R2/g) --gcBias -o salmon_GC/$(basename $f | cut - d. -f1 | sed s/_L001_R1//g) &&";done

library(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(DESeq2)
library(apeglm)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyr)
library(ggplot2)
library(fgsea)
library(tidyverse)
library(dplyr)

###Read in the counts from each sample listed in the sample table
samples <- read.table("DNMT1_AML_RNAseq_THP1_samples.txt", header = TRUE)
samples
files <- file.path("path_to_salmon_output", samples$run, "quant.sf")
names(files) <- paste0("sample", 1:52)
file.exists(files)

###Annotate the transcripts
###You need to make the TxDb from gff if the download is from Gencode
Txdb <- makeTxDbFromGFF("path_to/gencode.GRCh38v23.annotation.gtf", format="auto")
k <- keys(Txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(Txdb, k, "GENEID", "TXNAME")
head(tx2gene)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)

###Read the annotated counts table into DEseq
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Replicate + Treatment)

###Filter very low expressing across all samples genes 
keep <- rowSums(counts(ddsTxi)) >= 2
dds <- ddsTxi[keep,]

###Indicate which are the control samples
dds$experiment <- relevel(dds$Treatment, ref = "DMSO")

dds$group <- factor(paste0(dds$Dose,".", dds$Treatment,".",dds$Time))
dds$group
design(dds) <- ~ Replicate + group

###Perform differential expression
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

####Generate a TPM counts table###########
count<-counts(dds, normalized=TRUE)
colnames(count)<-t(samples$run)
genes<-data.frame(rownames(count))
gen<-genes %>% separate(rownames.count., c("Gene_ID", "database"))
rownames(count)<-gen$Gene_ID
write.table(count, file="DNMT1_AML_TPM_THP1_table.txt")

#Report results with annotation
res <-results(dds, contrast=c("group","400nM.DAC.hr144","0nM.DMSO.hr144"), alpha=0.05) 
###Replace the ESNG name source text for matching##############
genes<-data.frame(rownames(res))
gen<-genes %>% separate(rownames.res., c("Gene_ID", "database"))
rownames(res)<-gen$Gene_ID

#Get the stats from the differential
mcols(res, use.names=TRUE)
summary(res)

# Annotate the list of gene changes
res$Symbol <- mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
res$GeneName <- mapIds(org.Hs.eg.db,
                       keys = row.names(res),
                       column = "GENENAME",
                       keytype = "ENSEMBL",
                       multiVals = "first")

resOrdered <- res[order(res$padj),]
head(resOrdered)
# Print out results with annotations
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "DNMT1_THP1_400nMDAC_d6_DEseq2.csv")

####Perform GSEA using fGSEA######################
###Load the gene sets downloaded from MsigDB###
#gmt.file <- system.file("extdata", "h.all.v6.2.symbols.gmt", package="fgsea")
#gmt.file <- system.file("extdata", "custom_c2_h_all_symbols.gmt", package="fgsea")
gmt.file <- system.file("extdata", "Custom_genesets_hallmark_C2_symbols.gmt", package="fgsea")
pathways <- gmtPathways(gmt.file)
str(head(pathways))

###Get only the symbol and fold change from the differential result and get rid of NAs and average probes to for the same gene###
gene_list2 <- resOrderedDF %>% 
  dplyr::select(Symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Symbol) %>% 
  summarize(logFC=mean(stat))
gene_list2
ranks <- deframe(gene_list2)
head(ranks, 20)

fgseaRes <- fgsea(pathways = pathways, 
                  stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)

###sort the table###
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
pdf()
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
fgseaResTidy$gsea<-sapply(fgseaResTidy$leadingEdge,function(x) paste(x,collapse=","))
fgseaResTidy<-fgseaResTidy[,-8]
write.table(fgseaResTidy, file="DNMT1_THP1_400nM032_d4_fGSEA.txt")



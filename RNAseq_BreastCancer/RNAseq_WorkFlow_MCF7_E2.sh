Fastq -> Trimming -> Alignment -> Remove duplicate -> MACS2 peakcalling (--keepdup all) -> Turn to R


======================================================================================================================
#1. Download rawdata, QC and Triming   
workDir="/media/hkh/8TB/XUANTHANG/BreastCancer/MCF7_E2_RNAseq"
mkdir -p $workDir/RawData/sra $workDir/RawData/Fastq $workDir/RawData/FastQC $workDir/Results/STAR
export PATH=/media/hkh/8TB/XUANTHANG/TOOLS/sratoolkit.2.10.9-ubuntu64/bin:$PATH

RawData="/media/hkh/8TB/XUANTHANG/BreastCancer/MCF7_E2_RNAseq/RawData"
Results="/media/hkh/8TB/XUANTHANG/BreastCancer/MCF7_E2_RNAseq/Results"


# RNAseq SINGLE-END
prefetch --option-file $RawData/SRR*.txt --output-directory $RawData/sra/
fastq-dump $RawData/sra/SRR*/SRR* --outdir $RawData/Fastq --gzip
fastqc -t 6 -o $RawData/FastQC --noextract -f fastq $RawData/Fastq/SRR*
multiqc -o $RawData -n RNAseq_QC_report $RawData/FastQC/.

=====================================================================================================================
#2 Alignment RNAseq  (Single-end, STAR Using UCSC HUMAN hg38) 

# set up file names and grab base of filename for naming outputs

STARRef="/media/hkh/8TB/XUANTHANG/References/Refrence_Human_Genome/Homo_sapiens_UCSC_hg38/Sequence/STAR_Index/"
mkdir $Results/STARAlign

#Alignment multiple file (If many file, using --genomeLoad LoadAndKeep \)
for i in $RawData/Fastq/*.gz ; do
	echo $i;
	base=`basename -s .fastq.gz $i`
	STAR --runThreadN 7 \
	--genomeDir $STARRef \
	--runMode alignReads \
	--readFilesIn $RawData/Fastq/${base}.fastq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesCommand zcat \
	--outFileNamePrefix $Results/STARAlign/${base} 
done

# Quality check of aligned reads
multiqc $Results/STARAlign/. -o $Results/STARAlign -n STAR_alignment

#5  Downstream Analysis 
#5.1 RNAseq - Counting reads with SubRead
GeneAnnotation="/media/hkh/8TB/XUANTHANG/References/Refrence_Human_Genome/Homo_sapiens_UCSC_hg38/Annotation/Genes/genes.gtf"


mkdir $Results/feartureCounts
featureCounts \
      -t exon \
      -g gene_id \
      --primary \
      -a $GeneAnnotation \
      -o $Results/featureCounts $Results/STARAlign/*.bam 
#Then move to Rstudio
-----------------------








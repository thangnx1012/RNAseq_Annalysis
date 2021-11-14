#Example commands run for TETranscripts

### Using the unique alignments, generate a LTR family count table to input into DEseq2
bedtools multicov S1-09_S9_001.Aligned.out.bam S1-87_S87_001.Aligned.out.bam S1-10_S10_001.Aligned.out.bam S1-88_S88_001.Aligned.out.bam S1-11_S11_001.Aligned.out.bam S1-89_S89_001.Aligned.out.bam S1- 12_S12_001.Aligned.out.bam S1-90_S90_001.Aligned.out.bam S1-13_S13_001.Aligned.out.bam S1- 91_S91_001.Aligned.out.bam S1-14_S14_001.Aligned.out.bam S1-92_S92_001.Aligned.out.bam S1- 15_S15_001.Aligned.out.bam S1-93_S93_001.Aligned.out.bam -bed /media/epigenetics/01gsk/Annotations/RepeatMask_GRCh38_ERV.bed -s > DNMT1_MV411_d4_GSK032dose_counts.txt


### Align stranded total RNA-seq dataset to hg38 using --outFilterMultimapNmax ###100 for multimapping alignments to LTR retroelements###

for f in fastq/*L001_R1_001.fastq.gz; do echo "STAR --runThreadN 16 --genomeDir GRCh38/ --readFilesIn $f,$(echo $f | sed s/_L001_R1/_L002_R1/g),$(echo $f | sed s/_L001_R1/_L003_R1/g),$(echo $f | sed s/_L001_R1/_L004_R1/g),$(echo $f | sed s/_L001_R1/_L005_R1/g),$(echo $f | sed s/_L001_R1/_L006_R1/g) $(echo $f | sed s/L001_R1/L001_R2/g),$(echo $f | sed s/_L001_R1/_L002_R2/g),$(echo $f | sed s/_L001_R1/_L003_R2/g),$(echo $f | sed s/_L001_R1/_L004_R2/g),$(echo $f | sed s/_L001_R1/_L005_R2/g),$(echo $f | sed s/_L001_R1/_L006_R2/g) --sjdbGTFfile GRCh38/HS_gencodeGRCh38.v23.geneidedit.annotation.gtf --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix multimap/$(basename $f | cut -d. -f1 | sed s/_L001_R1//g). &&";done

###Remove the rRNA and mitochondria (chrM) from the alignments
split_bam.py -i S1-01_S1_001.Aligned.sortedByCoord.out -r RseQC/GRCh38_rRNA_chM.bed -o DNMT1_01_MV411_DMSO_d1_n1_noM

###Using the multimapping alignments, generate a count table and call differential ###expression using TEtoolkit###
TEtranscripts --format BAM --sortByPos --stranded reverse -t DNMT1_12_MV41 1_80nM032_d4_n1_sorted.ex.bam DNMT1_90_MV41 1_80nM032_d4_n2_sorted.ex.bam -c DNMT1_09_MV411_DMSO_d4_n1_sorted.ex.bam DNMT1_87_MV41 1_DMSO_d4_n2_sorted.ex.bam --GTF HS_gencodeGRCh38.v23.geneidedit.annotation.gtf --TE GRCh38_rmsk_TE.gtf --mode multi --project DNMT1_MV411_80nM032_d4_TE

### Separate the forward and reverse alignments for generating browser tracks###
samtools view -b -f 128 -F 16 input.bam > sample_name_fwd1.bam
samtools view -b -f 80 input.bam > sample_name_fwd2.bam
samtools merge -f sample_name_fwd.bam sample_name_fwd1.bam sample_name_fwd2.bam
samtools view -b -f 144 input.bam > sample_name_rev1.bam
samtools view -b -f 64 -F 16 input.bam > sample_name_rev2.bam
samtools merge -f sample_name_rev.bam sample_name_rev1.bam sample_name_rev2.bam

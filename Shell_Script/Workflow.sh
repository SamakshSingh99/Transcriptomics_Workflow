


# Generating Genome Index Files

STAR --runMode genomeGenerate \
--genomeDir ./Genome_Dir \
--genomeFastaFiles ./Reference/GCA_015475615.1_ASM1547561v1_genomic.fna \
--sjdbGTFfile ./Reference/GCA_015475615.1_ASM1547561v1_genomic.gtf \
--runThreadN 8

# Running Star Alignment 

STAR --genomeDir ./Genome_Dir \
--readFilesIn ./Dat_Dir/file_1.fastq.gz ./Dat_Dir/file_2.fastq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix ./results/wt2_ \
--outSAMtype BAM SortedByCoordinate \
--sjdbGTFfile ./Reference/GCA_015475615.1_ASM1547561v1_genomic.gtf \
--runThreadN 8

# Running HTSeq for gene counts

htseq-count -s no -f bam \
./results/file_Aligned.sortedByCoord.out.bam \
./Reference/GCA_015475615.1_ASM1547561v1_genomic.gtf \
> ./results/counts_file/file_count.txt
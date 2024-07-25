# RNA-Seq Analysis Pipeline: From Raw Reads to Gene Counts
This repository contains a command-line pipeline for analyzing RNA sequencing (RNA-Seq) data. The pipeline uses HISAT2 for aligning reads to a reference genome and HTSeq-count for quantifying gene expression levels.

## Pipeline Steps:

  ### 1. Create Directories: 
  >**mkdir -p sam_file bam_file count_file genome_index:** Sets up the necessary folders to organize            
            the output files at different stages of the pipeline.

  ### 2. Unzip FASTQ Files:
  >**gunzip ./raw_data/*.fastq.gz:** Extracts the raw sequencing data (FASTQ files) from their compressed 
            format (.fastq.gz) and places them in the raw_data directory. 

  ### 3. Build Genome Index:
  >**hisat2-build ./refgenome/hg38.fa ./genome_index/hg38_index:** HISAT2 is a fast and efficient aligner. This command creates an index of the reference genome (hg38 in this case) to speed up the alignment process. The index is stored in the genome_index directory.

  ### 4. Trim Reads (Quality Control):
  >**for file_R1 in ./raw_data/*_R1.fastq; do ... done:** This loop iterates over each paired-end FASTQ file in the raw_data directory.<br>
  >**trim_galore --paired "$file_R1" "$file_R2" --fastqc -o ./trimmed_fastq/:** Trim Galore is used for quality trimming (removing low-quality bases) and adapter removal from the reads. The --fastqc option generates FastQC reports for quality assessment. The trimmed reads are saved in the trimmed_fastq directory.

  ### 5. Align Reads (HISAT2):
  >**for file_R1 in ./trimmed_fastq/*_R1_val_1.fq; do ... done:** This loop iterates over the trimmed FASTQ files.<br>
  >**hisat2 --dta-cufflinks -x ./genome_index/hg38_index -1 "$file_R1" -2 "$file_R2" -S "$output_name":** HISAT2 aligns the trimmed reads to the reference genome using the previously built index. The --dta-cufflinks option prepares the alignments for downstream analysis with Cufflinks (a tool for transcript assembly and quantification). The aligned reads are saved in SAM format in the sam_file directory.

  ### 6. Convert to BAM:
  >**for converted_sam in ./sam_file/*.sam; do ... done:** This loop converts the SAM files (text-based alignment format) into BAM files (binary format) using samtools view. BAM files are smaller and more efficient for further analysis.

  ### 7. Count Reads (HTSeq-count):
  >**htseq-count -m union -f bam --additional-attr=transcript_id -s yes ./bam_file/*.bam ./gtf/hg38.ncbiRefSeq.gtf > count_file/htcount.txt:** HTSeq-count quantifies the number of reads that align to each gene or transcript in the genome. The -m union option specifies how to handle reads that overlap multiple features. The gene/transcript information is obtained from the GTF (Gene Transfer Format) file, which contains annotations about the genome.
>

  ### 8. Clean Up Count File:
  >**sed '/^__/ d' < count_file/htcount.txt > count_file/final_htcount.txt:** Removes lines starting with "__" from the HTSeq-count output. These lines typically contain summary information, not individual gene counts. The cleaned count data is saved as final_htcount.txt.<br>

## Usage
1. Clone this repository (git clone [repository URL])
2. Install dependencies
3. Place your raw FASTQ files in the raw_data directory.
4. Place your reference genome file (e.g., hg38.fa) in the refgenome directory.
5. Place your GTF annotation file (e.g., hg38.ncbiRefSeq.gtf) in the gtf directory.
6. Execute the pipeline script: bash pipeline.sh (or the appropriate command for your shell)

## Dependencies
<pre>
* <a href="https://anaconda.org/bioconda/hisat2">HISAT2</a>
* <a href="https://anaconda.org/bioconda/samtools">Samtools</a>
* <a href="https://anaconda.org/bioconda/htseq">HTSeq</a>
* <a href="https://anaconda.org/bioconda/trim-galore">Trim Galore!</a>
* <a href="https://anaconda.org/bioconda/fastqc">FastQC</a>
</pre>

## Copyright
© [2024] Mohammad Uzzal Hossain. All rights reserved.

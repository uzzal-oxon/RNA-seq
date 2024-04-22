#!/bin/bash

# Clear and informative header
echo "-----------------------------------------"
echo "  RNA-Seq Analysis Pipeline Script"
echo "-----------------------------------------"

# Function to calculate time difference
calculate_time() {
    end=$1
    start=$2
    diff=$((end-start))
    echo "Time taken: $diff seconds"
}

# Create directories with error checking
start=$(date +%s)
for dir in genome_index trimmed_fastq SAM BAM features; do
    if ! mkdir "$dir"; then
        echo "Error: Failed to create directory $dir. Exiting."
        exit 1
    fi
done
end=$(date +%s)
calculate_time $end $start

# Build genome index
echo "Building genome index..."
start=$(date +%s)
hisat2-build -p 90 ./refgenome/hg38.fa genome_index/genome_index || {
    echo "Error: HISAT2 index build failed. Exiting."
    exit 1
}
end=$(date +%s)
calculate_time $end $start

# Trim FASTQ files
echo "Trimming FASTQ files..."
start=$(date +%s)
for file_R1 in ./sra/*_1.fastq; do
    file_R2="${file_R1%_1.fastq}_2.fastq"
    trim_galore --paired "$file_R1" "$file_R2" --fastqc -o ./trimmed_fastq/
done
end=$(date +%s)
calculate_time $end $start

# Align reads 
echo "Aligning reads to genome (This may take a while)..."
start=$(date +%s)
for file_R1 in ./trimmed_fastq/*_1_val_1.fq; do
    file_R2="${file_R1%_1_val_1.fq}_2_val_2.fq"
    output_name="./SAM/$(basename "${file_R1%_1_val_1.fq}").sam"
    hisat2 -p 180 --dta-cufflinks -x ./genome_index/genome_index -1 "$file_R1" -2 "$file_R2" -S "$output_name"
done
end=$(date +%s)
calculate_time $end $start

# Convert SAM to BAM 
echo "Converting SAM to BAM..."
start=$(date +%s)
for sam_file in ./SAM/*.sam; do
    output_bam="./BAM/$(basename "${sam_file%.sam}.bam")"
    samtools view -bS "$sam_file" > "$output_bam"
done
end=$(date +%s)
calculate_time $end $start

# Feature counting 
echo "Performing feature counting..."
start=$(date +%s)
htseq-count -m union -f bam --additional-attr=transcript_id -s yes ./BAM/*.bam ./gtf/hg38.ncbiRefSeq.gtf > features/htcount.txt

# Filter count results
echo "Filtering count results..."
sed '/^__/ d' < features/htcount.txt > features/htcount2.txt

end=$(date +%s)
calculate_time $end $start

echo "-----------------------------------------"
echo "Pipeline Complete!"
echo "-----------------------------------------"


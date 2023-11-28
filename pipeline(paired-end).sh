#!/bin/bash

# Specify the main directory where datasets are located
MAIN_DIR="$HOME/Desktop/Prakhar"

# Specify the names of your datasets
DATASETS=("Mention the dataset name here")
SRA=("Mention the SRA numbers here")

# Function to print elapsed time
print_elapsed_time() {
    local start_time=$1
    local end_time=$2
    echo "Time taken: $((end_time - start_time)) seconds"
}

# Iterate over datasets
for dataset in "${DATASETS[@]}"; do
    INPUT_DIR="$MAIN_DIR/data/$dataset"
    FASTQC_DIR="$MAIN_DIR/fastq-results/$dataset"
    TRIM_DIR="$MAIN_DIR/trim/$dataset"
    HISAT_DIR="$MAIN_DIR/hisat-results/$dataset"
    GENOME_DIR="$MAIN_DIR/genome"
    SAMTOOLS_DIR="$HISAT_DIR"  # Samtools output will be in the same directory as HISAT2
    STRINGTIE_DIR="$MAIN_DIR/stringtie/$dataset"

    # Create output directories if they don't exist
    mkdir -p "$INPUT_DIR"
    mkdir -p "$FASTQC_DIR"
    mkdir -p "$TRIM_DIR"
    mkdir -p "$HISAT_DIR"
    mkdir -p "$STRINGTIE_DIR"

    # Iterate over sample files in the dataset
    for file_1 in "${SRA[@]}"; do
        # Check if the file is a regular file

        sample_name=$file_1

        #Downloading the files
        echo "----------Downloading $sample_name in $dataset----------"
        start_time_fastqc=$(date +%s)
        prefetch -p "$sample_name" -O "$INPUT_DIR" && fasterq-dump -pe 24 "$INPUT_DIR/$sample_name/$sample_name".sra -O "$INPUT_DIR"
        end_time_fastqc=$(date +%s)
        print_elapsed_time "$start_time_fastqc" "$end_time_fastqc"

        #Deleting the sra files
        rm -r "$INPUT_DIR/$sample_name"

        # Run FastQC
        echo "----------Running FastQC for $sample_name in $dataset----------"
        start_time_fastqc=$(date +%s)
        fastqc "$INPUT_DIR/$sample_name"_1.fastq "$INPUT_DIR/$sample_name"_2.fastq -o "$FASTQC_DIR"
        end_time_fastqc=$(date +%s)
        print_elapsed_time "$start_time_fastqc" "$end_time_fastqc"

        # Run Trimmomatic
        echo "----------Running Trimmomatic for $sample_name in $dataset----------"
        start_time_trimmomatic=$(date +%s)
        mkdir -p "$TRIM_DIR/$sample_name"
        java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 "$INPUT_DIR/$sample_name"_1.fastq "$INPUT_DIR/$sample_name"_2.fastq.gz -baseout "$TRIM_DIR/$sample_name/$sample_name.fastq.gz" ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE.fa:4:30:10 MINLEN:36 TRAILING:3 LEADING:3
        end_time_trimmomatic=$(date +%s)
        print_elapsed_time "$start_time_trimmomatic" "$end_time_trimmomatic"

        #Deleting the fastq files
        rm "$INPUT_DIR/$sample_name".fastq

        # Run hisat2
        echo "----------Running hisat2 for $sample_name in $dataset----------"
        start_time_hisat2=$(date +%s)
        mkdir -p "$HISAT_DIR/$sample_name"
        hisat2 -q -x "$GENOME_DIR/genome" -1 "$TRIM_DIR/$sample_name/$sample_name"_1P.fastq.gz -2 "$TRIM_DIR/$sample_name/$sample_name"_2P.fastq.gz -S "$HISAT_DIR/$sample_name/$sample_name.sam" -p 16
        end_time_hisat2=$(date +%s)
        print_elapsed_time "$start_time_hisat2" "$end_time_hisat2"

        # Delete trimmed files after hisat2
        rm "$TRIM_DIR/$sample_name/$sample_name"*fastq.gz
        

        # Run samtools
        echo "----------Running samtools for $sample_name in $dataset----------"
        start_time_samtools=$(date +%s)
        samtools view -bS "$HISAT_DIR/$sample_name/$sample_name.sam" -o "$HISAT_DIR/$sample_name/$sample_name.bam"
        samtools sort "$HISAT_DIR/$sample_name/$sample_name.bam" -o "$HISAT_DIR/$sample_name/$sample_name.sorted.bam"
        mv "$HISAT_DIR/$sample_name/$sample_name.sorted.bam" "$HISAT_DIR/$sample_name/final_merged.bam"
        samtools index "$HISAT_DIR/$sample_name/final_merged.bam"
        end_time_samtools=$(date +%s)
        print_elapsed_time "$start_time_samtools" "$end_time_samtools"

        # Delete SAM files after samtools
        rm "$HISAT_DIR/$sample_name/$sample_name.sam"

        # Run StringTie 
        echo "----------Running StringTie for $sample_name in $dataset----------"
        start_time_stringtie=$(date +%s)
        mkdir -p "$STRINGTIE_DIR/$sample_name"
        stringtie "$HISAT_DIR/$sample_name/final_merged.bam" -G "$GENOME_DIR/genome.gtf" -e -o "$STRINGTIE_DIR/$sample_name/$sample_name.gtf" -p 24 -l "$sample_name"
        end_time_stringtie=$(date +%s)
        print_elapsed_time "$start_time_stringtie" "$end_time_stringtie"

        #Delete the bam files folder
        rm -r "$HISAT_DIR/$sample_name"

        # Create a text file with Sample_name and sample_gtf_file_location
        echo "$sample_name $STRINGTIE_DIR/$sample_name/$sample_name.gtf" >> "$STRINGTIE_DIR/sample_info.txt"

        # Print time taken for the entire analysis for this sample
        echo "----------Total Time taken for $sample_name in $dataset----------"
        print_elapsed_time "$start_time_fastqc" "$end_time_stringtie"
    done
done


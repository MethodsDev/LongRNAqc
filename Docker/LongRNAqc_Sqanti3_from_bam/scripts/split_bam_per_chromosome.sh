#!/usr/bin/env bash

# split BAM per chromosome

if [ $# -ne 4 ]; then
    echo "Usage: $0 input_bam output_folder bam|sam chromosomes"
    exit 1
fi

input_bam="$1"
output_folder="$2"
file_type=$(echo "$3" | tr '[:upper:]' '[:lower:]')  # Convert to lowercase for case insensitivity
chromosomes_input="$4"

if [[ "$file_type" != "bam" && "$file_type" != "sam" ]]; then
    echo "Invalid file type. Must be either 'bam' or 'sam'."
    exit 1
fi

if [ ! -f "$input_bam" ]; then
    echo "Input BAM file not found: $input_bam"
    exit 1
fi

if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

# Extract chromosome names from the BAM file
bam_chromosomes=$(samtools view -H "$input_bam" | grep -E '^@SQ' | cut -f 2 | cut -d ':' -f 2 | sort -u)

# Convert the comma-separated input into an array and then intersect with bam chromosomes
IFS=',' read -ra chrom_array <<< "$chromosomes_input"
intersected_chromosomes=()
for chrom in "${chrom_array[@]}"; do
    if echo "$bam_chromosomes" | grep -qw "$chrom"; then
        intersected_chromosomes+=("$chrom")
    fi
done

# Split per chromosome only for the intersection
for chromosome in "${intersected_chromosomes[@]}"; do
    if [ "$file_type" == "bam" ]; then
        output_file="${output_folder}/${chromosome}.bam"
        samtools view -b "$input_bam" "$chromosome" > "$output_file"
    else
        output_file="${output_folder}/${chromosome}.sam"
        samtools view "$input_bam" "$chromosome" > "$output_file"
    fi
done

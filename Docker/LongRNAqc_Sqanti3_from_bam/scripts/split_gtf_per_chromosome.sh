#!/usr/bin/env bash

# split GTF per chromosome

if [ $# -ne 3 ]; then
    echo "Usage: $0 input_gtf output_folder chromosomes"
    exit 1
fi

input_gtf="$1"
output_folder="$2"
chromosomes_input="$3"

if [ ! -f "$input_gtf" ]; then
    echo "Input GTF file not found: $input_gtf"
    exit 1
fi

if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

# Extract chromosome names from the GTF file
bam_chromosomes=$(grep -v '^#' "$input_gtf" | cut -f1 | sort -u)

# Convert the comma-separated input into an array and then intersect with chromosomes from GTF
IFS=',' read -ra chrom_array <<< "$chromosomes_input"
intersected_chromosomes=()
for chrom in "${chrom_array[@]}"; do
    if echo "$bam_chromosomes" | grep -qw "$chrom"; then
        intersected_chromosomes+=("$chrom")
    fi
done

# Split per chromosome only for the intersection
for chromosome in "${intersected_chromosomes[@]}"; do
    output_file="${output_folder}/${chromosome}.gtf"
    awk -v chrom="$chromosome" '!/^#/ && $1 == chrom' "$input_gtf" > "$output_file"
done

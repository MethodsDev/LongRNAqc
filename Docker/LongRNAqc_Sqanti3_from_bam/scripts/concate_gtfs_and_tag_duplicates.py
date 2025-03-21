#!/usr/bin/env python3

import argparse
import re

# usage: $ python concate_gtfs_and_tag_duplicates.py -o combined_output.gtf in_file1.gtf in_file2.gtf

def process_files(input_files, output_file):
    pattern = re.compile(r'_dup\d+$')
    id_counts = {}
    with open(output_file, 'w') as out_f:
        for file in input_files:
            with open(file, 'r') as in_f:
                for line in in_f:
                    cols = line.strip().split('\t')
                    if cols[2] == "transcript":
                        transcript_id = (cols[8].split(";")[0])[15:-1] # remove starting 'transcript_id "' and trailing '"'
                        transcript_id = pattern.sub('', transcript_id) # remove any potential "_dupN" that already exists
                        if transcript_id in id_counts:
                            id_counts[transcript_id] += 1
                            current_id = f"{transcript_id}_dup{id_counts[transcript_id]}"
                        else:
                            id_counts[transcript_id] = 1
                            current_id = transcript_id
                        cols[8] = f'transcript_id "{current_id}"; gene "{current_id}"'
                    elif cols[2] == "exon":
                        cols[8] = f'transcript_id "{current_id}"'
                    out_f.write('\t'.join(cols) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GTF files and combine them, appending _dupN to multimapping reads.\n\
    usage: $ python concate_gtfs_and_tag_duplicates.py -o combined_output.gtf in_file1.gtf in_file2.gtf")
    parser.add_argument('input_files', metavar='FILES', type=str, nargs='+',
                        help='List of input GTF files to be processed')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Output GTF file where the results will be written')

    args = parser.parse_args()
    process_files(args.input_files, args.output)

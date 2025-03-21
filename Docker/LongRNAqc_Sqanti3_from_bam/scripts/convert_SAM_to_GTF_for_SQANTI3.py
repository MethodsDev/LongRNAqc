#!/usr/bin/env python3

import os, sys, subprocess
from Bio import SeqIO
from err_correct_w_genome import err_correct
from sam_to_gff3 import convert_sam_to_gff3
import argparse

# GFFREAD_PROG = "gffread"

def main():
    parser = argparse.ArgumentParser(description="Convert SAM input to a GTF file, optionally generating a corrected FASTA file of reads as well.")
    parser.add_argument("--reference_genome", required=True, help="Path to the reference genome FASTA file")
    parser.add_argument("--sam_file", required=True, help="Path to the corrected SAM file")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output corrected FASTA and GTF files")
    parser.add_argument("--allow_non_primary", action="store_true", help="Process non primary alignments as well")
    parser.add_argument("--fasta_correction", action="store_true", help="Perform FASTA correction using err_correct")
    args = parser.parse_args()


    corrGTF = args.output_prefix + ".gtf"

    if args.fasta_correction:
        corrFASTA = args.output_prefix + ".corrected.fasta"
        genome_dict = dict((r.name, r) for r in SeqIO.parse(open(args.reference_genome), 'fasta'))
        err_correct(args.reference_genome, args.sam_file, corrFASTA, genome_dict=genome_dict)

    # in corrSAM, out corrGTF
    convert_sam_to_gff3(args.sam_file, corrGTF, source=os.path.basename(args.reference_genome).split('.')[0], allow_non_primary=args.allow_non_primary)
    # cmd = "{p} {o}.tmp -T -o {o}".format(o=corrGTF, p=GFFREAD_PROG)
    # if subprocess.check_call(cmd, shell=True) != 0:
    #     print("ERROR running cmd: {0}".format(cmd), file=sys.stderr)
    #     sys.exit(-1)

if __name__ == "__main__":
    main()

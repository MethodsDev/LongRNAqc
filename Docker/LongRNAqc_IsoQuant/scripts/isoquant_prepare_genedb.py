#!/usr/bin/env python3

import argparse
import json
import os.path
import sys
from src.gtf2db import convert_gtf_to_db

def main(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--clean_start', action='store_true', default=False,
                        help='Do not use previously generated index, feature db or alignments.')
    parser.add_argument('--complete_genedb', action='store_true', default=False,
                        help="use this flag if gene annotation contains transcript and gene metafeatures, "
                             "e.g. with official annotations, such as GENCODE; "
                             "speeds up gene database conversion")
    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF format", type=str)
    parser.add_argument("--genedb_output", help="output folder for converted gene database,"
                                                  " will be created automatically (same as output by default)", type=str)
    args = parser.parse_args(args)

    config_dir = os.path.join(os.environ['HOME'], '.config', 'IsoQuant')
    os.makedirs(config_dir, exist_ok=True)

    args.gtf_check = True
    args.db_config_path = os.path.join(config_dir, 'db_config.json')
    args.index_config_path = os.path.join(config_dir, 'index_config.json')
    args.bed_config_path = os.path.join(config_dir, 'bed_config.json')
    args.alignment_config_path = os.path.join(config_dir, 'alignment_config.json')
    for config_path in (args.db_config_path, args.index_config_path, args.bed_config_path, args.alignment_config_path):
        if not os.path.exists(config_path):
            with open(config_path, 'w') as f_out:
                json.dump({}, f_out)

    convert_gtf_to_db(args)

if __name__ == "__main__":
    main(sys.argv[1:])

#!/usr/bin/env python

import argparse
from collections import defaultdict

intronic_counts = {}
all_counts = {}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to deduplicate sciRNA data given stdin of sorted BAM file. BAM printed to STDOUT. Tolerates 1bp mismatches in UMIs.')
    parser.add_argument('--gene_assignment_files', required=True, nargs='+', help='List of input files with gene assignments.') 
    parser.add_argument('--all_counts_file', required=True, help='Counts per cell in exons.') 
    parser.add_argument('--intron_counts_file', required=True, help='Counts per cell in introns.') 
    args = parser.parse_args()

    for f in args.gene_assignment_files:
        for line in open(f):
            entries = line.strip().split('\t')
            
            category = entries[2]
            barcodes = entries[0].split('|')
            sample, cell = barcodes[0], barcodes[1]

            if category == 'exonic' or category == 'intronic':
                cell_key = (sample, cell)

                all_counts[cell_key] = all_counts.get(cell_key, 0) + 1

                if category == 'intronic':
                    intronic_counts[cell_key] = intronic_counts.get(cell_key, 0) + 1
                    
    with open(args.all_counts_file, 'w') as all_counts_file, open(args.intron_counts_file, 'w') as intron_counts_file:
        for count_key in all_counts:
            all_count = all_counts.get(count_key, 0)
            intron_count = intronic_counts.get(count_key, 0)
            
            sample, cell = count_key
            all_counts_file.write(f'{sample}\t{cell}\t{all_count}\n')
            intron_counts_file.write(f'{sample}\t{cell}\t{intron_count}\n')

#!/usr/bin/env python

import argparse
from collections import defaultdict

intronic_counts = {}
all_counts = {}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to count umis per cell.')
    parser.add_argument('--gene_assignment_files', required=True, nargs='+', help='List of input files with gene assignments.') 
    parser.add_argument('--all_counts_file', required=True, help='Counts per cell in exons and introns.') 
    parser.add_argument('--intron_counts_file', required=True, help='Counts per cell in introns.') 
    args = parser.parse_args()

    for f in args.gene_assignment_files:
        for line in open(f):
            entries = line.strip().split('\t')
            
            category = entries[2]
            cell = entries[0]

            if category == 'exonic' or category == 'intronic':
                all_counts[cell] = all_counts.get(cell, 0) + 1

                if category == 'intronic':
                    intronic_counts[cell] = intronic_counts.get(cell, 0) + 1
                    
    with open(args.all_counts_file, 'w') as all_counts_file, open(args.intron_counts_file, 'w') as intron_counts_file:
        for cell in all_counts:
            all_count = all_counts.get(cell, 0)
            intron_count = intronic_counts.get(cell, 0)
            
            all_counts_file.write(f'{cell}\t{all_count}\n')
            intron_counts_file.write(f'{cell}\t{intron_count}\n')

#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import pysam
from copy import copy

BASES = ['A', 'T', 'G', 'C', 'N']

def generate_1bp_mismatches(sequence):
    mismatches = []
    sequence = list(sequence)

    for i in range(len(sequence)):
        for base in BASES:
            if base != sequence[i]:
                new_string = copy(sequence)
                new_string[i] = base
                mismatches.append(''.join(new_string))

    return mismatches

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to deduplicate sciRNA data given sorted BAM file. BAM printed to STDOUT. Tolerates 1bp mismatches in UMIs.')
    parser.add_argument('--bam', help='Input BAM file to be deduplicated.', required=True) 
    parser.add_argument('--output_bam', help='Output BAM file, deduplicated.', required=True)
    args = parser.parse_args()

    cell_umis_at_position = set() # (cell, umi)
    current_position = None

    readsin = pysam.AlignmentFile(args.bam, "rb")
    readsout = pysam.AlignmentFile(args.output_bam, "wb", template=readsin)
    reads_included = 0
    position_reads = dict()
    for read in readsin:
        position = read.reference_start
        if not current_position or current_position == position:
            position_reads[read.query_name] = read
            reads_included += 1
            current_position = position
        else:
            for key,value in sorted(position_reads.items()):
                read_name_parts = key.split('|')
                cell_barcode = f'{read_name_parts[2]}_{read_name_parts[3]}_{read_name_parts[4]}'
                umi = read_name_parts[5]
                cell_umi_key = (cell_barcode, umi)
       	        if cell_umi_key not in cell_umis_at_position:
                    readsout.write(value)               
                # Track the UMI and all 1bp mismatches to it
                cell_umis_at_position.add((cell_barcode, umi))

                for mismatch in generate_1bp_mismatches(umi):
                    cell_umis_at_position.add((cell_barcode, mismatch))
            cell_umis_at_position = set()
            position_reads = dict()
            position_reads[read.query_name] = read
            reads_included += 1
            current_position = position
    print(reads_included)

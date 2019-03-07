#!/usr/bin/env python
# Make sample fastqs 2-level
# Andrew's barcode parser/fastq combiner modified for 2-level

import barcodeutils as bu
import argparse
import os
import json
from collections import OrderedDict
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
P5_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p5.txt')
P7_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p7.txt')
RT_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/rt2.txt')



def get_programmed_pcr_combos(p5_lookup, p7_lookup, p5_cols_used, p7_rows_used):
    """
    Assuming p5 are columns and p7 are rows, get the set of (p5, p7) wells that were programmed.
    These are the set that would have been physically combined.
    Args:
        p5_lookup (dict): p5_lookup dict mapping sequences to wells as passed to barcode specification
        p7_lookup (dict): p7_lookup dict mapping sequences to wells as passed to barcode specification
        p5_cols_used (list): A list of the cols used from P5 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. A B C for p7 and 1 2 3 for p5.
        p7_rows_used: A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. A B C for p7 and 1 2 3 for p5.
    Returns:
        set of (p5, p7): set of (p5, p7) well ID tuples that are valid.
    """

    valid_combos = set()

    for p5, p7 in zip(p5_cols_used, p7_rows_used):

        selected_p5 = [p5_well for p5_well in p5_lookup.values() if int(p5_well[1:]) == p5]
        selected_p7 = [p7_well for p7_well in p7_lookup.values() if p7_well[0] == p7]

        for selected_p5_i in selected_p5:
            for selected_p7_i in selected_p7:
                valid_combos.add((selected_p5_i, selected_p7_i))

    return valid_combos


def quick_parse(file_path):
    # a copy of only the relevant lines from easygrid.read_delim
    fh = open(file_path)
    columns = next(fh).strip().split(",")

    # Parse file
    for line_number, line in enumerate(fh):
        entries = line.strip().split(",")

        if len(entries) != len(columns):
            raise ValueError('Length of entries: %s, not equal to length of columns %s, in file: %s, line number %s' % (len(entries), len(columns), file_path, line_number))

        entries_dict = dict(zip(columns, entries))
        yield entries_dict

def load_sample_layout(file_path):
    """
    Function that loads the sample layout file to an RT p5_lookup table.
    """

    lookup = {}
    for rt_well in quick_parse(file_path):
        lookup[rt_well['RT Barcode']] = rt_well['Sample ID'].replace('-', '.').replace('_', '.').replace(' ', '.')
    return lookup 




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to combine R1 and R2 file from sci run into a single file with barcodes in the fastq header.')
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    parser.add_argument('--read1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
    parser.add_argument('--read2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
    parser.add_argument('--file_name', required=True, help='The R1 file name.')
    parser.add_argument('--sample_layout', required=True, help='Text file containing the sample layout by RT well.')
    parser.add_argument('--p5_cols_used', nargs='+', type=int, required=True, help='A list of the columns used from P5 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p7_rows_used', nargs='+', required=True, help='A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--output_dir', required=True, help='Output directory for files.')
    parser.add_argument('--p7_length', type=int, default=10, help='Expected P7 index length.')
    parser.add_argument('--p5_length', type=int, default=10, help='Expected P5 index length.')

    args = parser.parse_args()

    lane_num = args.file_name
    lane_num = lane_num.replace("Undetermined_S0_L", "L")
    lane_num = lane_num.replace("_R1_001.fastq.gz", "")
    stats_file = os.path.join(args.output_dir, lane_num + ".stats.json")
    suffix = lane_num + ".fastq"

    reverse_complement_i5 = bu.reverse_complement_i5(args.run_directory)

    run_info = bu.get_run_info(args.run_directory)

    # Load barcodes
    p7_lookup = bu.load_whitelist(P7_FILE)
    p7_lookup = {sequence[0:args.p7_length]: well for sequence,well in p7_lookup.items()}

    p5_lookup = bu.load_whitelist(P5_FILE)
    if reverse_complement_i5:
        p5_lookup = {bu.reverse_complement(sequence): well for sequence,well in p5_lookup.items()}
    p5_lookup = {sequence[0:args.p5_length]: well for sequence,well in p5_lookup.items()}

    # Get the set of all valid PCR combos
    programmed_pcr_combos = get_programmed_pcr_combos(p5_lookup, p7_lookup, args.p5_cols_used, args.p7_rows_used)

    # Define where all sequences are and what the whitelists are
    barcode_spec = {
        'p5': {
            'start': 1,
            'end':  args.p5_length,
            'read': 'i5',
            'whitelist': p5_lookup
        },
        'p7': {
            'start': 1,
            'end': args.p7_length,
            'read': 'i7',
            'whitelist': p7_lookup
        },
        'umi': {
            'start': 1,
            'end': 8,
            'read': 'r1'
        },
        'rt': {
            'start': 9,
            'end': 18,
            'read': 'r1',
            'whitelist': RT_FILE
        }  
    }


  # Set up the output files
    sample_rt_lookup = load_sample_layout(args.sample_layout)
    
    sample_to_output_filename_lookup = {sample: os.path.join(args.output_dir, '%s-%s' % (sample, suffix)) for well,sample in sample_rt_lookup.items()}
    sample_to_output_file_lookup = {sample: open(filename, 'w') for sample,filename in sample_to_output_filename_lookup.items()}
    print("Demuxing %s samples (%s total RT wells) into their own files..." % (len(sample_to_output_filename_lookup), len(sample_rt_lookup)))
     
    # Set up some basic tracking for each sample
    sample_read_counts = {}
    for sample in sample_to_output_file_lookup:
        sample_read_counts[sample] = 0

    total_reads = 0
    total_uncorrected = 0
    total_pcr_mismatch = 0
    total_unused_rt_well = 0
    total_corrected = 0

    # Finally, process reads
    for read_number, entry in enumerate(bu.parse_fastq_barcodes(args.read1, args.read2, spec=barcode_spec, edit_distance=1)):

        total_reads += 1

        # Only allow the programmed PCR combos (helps clean things up a bit)
        p5 = entry['p5']
        p7 = entry['p7']
 
        corrected = entry['rt'] is not None
        corrected_p5_p7 = p5 is not None and p7 is not None

        if corrected and corrected_p5_p7:
            total_corrected += 1
            rt_barcode = entry['rt']
            umi = entry['umi']
        else:
            total_uncorrected += 1
            continue

        if not (p5, p7) in programmed_pcr_combos:
            total_pcr_mismatch += 1
            continue

        if rt_barcode not in sample_rt_lookup:
            total_unused_rt_well += 1
            continue

        sample = sample_rt_lookup[rt_barcode]
        sample_read_number = sample_read_counts[sample] + 1
        sample_read_counts[sample] += 1

        r2_qual = entry['r2_qual']
        r2_seq = entry['r2_seq']
        output_name = f'@{sample}-P7{p5}-P5{p7}_{sample_read_number}|{sample}|{p5}|{p7}|{rt_barcode}|{umi}'
        output_line = f'{output_name}\n{r2_seq}\n+\n{r2_qual}\n'
        sample_to_output_file_lookup[sample].write(output_line)

    # Close output files
    for f in sample_to_output_file_lookup.values():
        f.close()

    # Output stats
    total_passed_reads = sum(list(sample_read_counts.values()))
    stats = OrderedDict()
    stats['total_input_reads'] = total_reads
    stats['total_passed_reads'] = total_passed_reads
    stats['fraction_passed_reads'] = total_passed_reads / total_reads
    stats['fraction_uncorrected_reads'] = total_uncorrected / total_reads
    stats['fraction_invalid_rt_well'] = total_unused_rt_well / total_reads
    stats['fraction_pcr_mismatch'] = total_pcr_mismatch / total_reads
    stats['total_reads_corrected'] = total_corrected
    stats['total_reads_passed_per_sample'] = sample_read_counts


    with open(stats_file, 'w') as f:
        f.write(json.dumps(stats, indent=4))


#!/usr/bin/env python
# Andrew's sample sheet creator

import os
import barcodeutils as bu
import argparse

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
P5_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p5.txt')
P7_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p7.txt')

def get_sample_sheet_text(p5_index_length, p7_index_length):
    """
    Gets the sample sheet text that will demux cells into one set of files
    """

    sample_sheet_template = """[DATA]
Lane,Sample_ID,Sample_Name,index,index2
%s"""

    line = ',fake,fake,%s,%s' % ('N' * p7_index_length, 'N' * p5_index_length)
    return sample_sheet_template % line


def load_barcode_file(barcode_file, index_length=None):
    """
    Loads one of the sciRNA barcode files allowing it to be trimmed to a specified length.
    """
    barcode_dict = {}
    for line in open(barcode_file):
        index, barcode = line.strip().split('\t')

        if index_length:
            if len(barcode) < index_length:
                raise ValueError('The specified barcodes are not long enough to match the run index length: %s vs. %s' % (len(barcode), index_length))
            barcode = barcode[0:index_length]

        if index in barcode_dict and barcode_dict[index] != barcode:
            raise ValueError('Conflicting entries in barcode file for index: %s, %s' % (barcode_file, index))
        barcode_dict[index] = barcode
    return barcode_dict


if __name__ == '__main__':
    # Script
    parser = argparse.ArgumentParser('Script make sample sheet')

    # Required args
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    args = parser.parse_args()

    # Get simple things like index lengths and flow cell ID for the run
    run_info = bu.get_run_info(args.run_directory)

    # Set up samplesheet for BCL2FASTQ
    p5_indices = load_barcode_file(P5_FILE)
    reverse_i5 = bu.reverse_complement_i5(args.run_directory)
    if reverse_i5:
        p5_indices = {x:bu.reverse_complement(y) for x,y in p5_indices.items()}
    p5_indices = {x: y[0:run_info['p5_index_length']] for x,y in p5_indices.items()}

    p7_indices = load_barcode_file(P7_FILE, run_info['p7_index_length'])

    sample_sheet_text = get_sample_sheet_text(run_info['p7_index_length'], run_info['p5_index_length'])
    sample_sheet_path = os.path.join('SampleSheet.csv')    
    with open(sample_sheet_path, 'w') as sample_sheet:
        sample_sheet.write(sample_sheet_text)


#!/usr/bin/env python    

import argparse
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
RT_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/rt2.txt')
RT3_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/rt.txt')

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Check sample sheet')

    parser.add_argument('--sample_sheet', required=True, help='Path to sample sheet.')
    parser.add_argument('--star_file', required=True, help='Path to star genomes file.')
    parser.add_argument('--level', required=True, help='2 or 3 level sci?')
    args = parser.parse_args()

    if args.level == "3":
        rtfile = RT3_FILE
    else:
        rtfile = RT_FILE

    rtdict = {}
    with open(rtfile) as rt_file:
        for line in enumerate(rt_file):
            line = line[1].strip().split("\t")
            rtdict[line[0]] = line[1]

    genomes = []

    with open(args.star_file, 'r') as f:
        for line in f:
            items = line.strip().split()
            genomes.append(items[0])
           
    def check_line(line, line_num, rtdict = rtdict, genomes = genomes):
        error_flag = 0
        line = line.strip().split(",")
        if not line[0] in rtdict.keys():
            sys.stderr.write("Sample sheet error at line " + str(line_num) + ". RT Barcode '" + line[0] + "' not valid.\n")
            error_flag = 1
        if not line[2] in genomes:
            sys.stderr.write("Sample sheet error at line " + str(line_num) + ". Reference Genome '" + line[2] + "' not valid.\n")
            error_flag = 1
        return error_flag
    sheet = open(args.sample_sheet)

    # check for header
    topline_orig = sheet.readline()
    topline = topline_orig.strip().split(",")
    line_num = 1
    num_cols = len(topline)
    if num_cols == 3:
        if topline[0] == 'RT Barcode' and topline[1] =='Sample ID' and topline[2] == 'Reference Genome':
            sample_out = 'RT Barcode,Sample ID,Reference Genome\n'
        else:
            sample_out = 'RT Barcode,Sample ID,Reference Genome\n'
            check_line(topline_orig, line_num)
            sample_out = sample_out + topline_orig
    else:
        if topline[0] == 'RT Barcode' and topline[1] =='Sample ID' and topline[2] == 'Reference Genome':
            sample_out = topline_orig + "\n"
        else:
            sample_out = 'RT Barcode,Sample ID,Reference Genome'
            for i in range(3, len(topline)):
                sample_out = sample_out + ",Column" + str(i)
            check_line(topline_orig, line_num)
            sample_out = sample_out + "\n" +  topline_orig
    # Check RT and Genomes against possible
    error_count = 0

    line_num = 1
    for line in sheet:
        line_num += 1
        linesp = line.strip().split(",")
        if linesp[0] in (None, "") and linesp[1] in (None, "") and linesp[2] in (None, ""):
            continue
        error_count += check_line(line, line_num)
        sample_out = sample_out + line
    sheet.close()
    
    if error_count > 0:
        sys.stderr.write("There were " + str(error_count) + " errors in the sample sheet.\n")
        sys.exit(10)

    new_sheet = open("good_sample_sheet.csv", "w")
    new_sheet.write(sample_out)
    new_sheet.close()

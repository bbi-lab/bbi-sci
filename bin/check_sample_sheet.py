#!/usr/bin/env python    

import argparse
import sys
import os
import math

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
RT_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/rt2.txt')
RT3_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/rt.txt')

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Check sample sheet')
    parser.add_argument('--sample_sheet', required=True, help='Path to sample sheet.')
    parser.add_argument('--star_file', required=True, help='Path to star genomes file.')
    parser.add_argument('--level', required=True, help='2 or 3 level sci?')
    parser.add_argument('--rt_barcode_file', required=True, help='Custom barcode file path or "default"')
    parser.add_argument('--max_wells_per_sample', required=True, help='Maximum number of wells per sample - for efficiency')
    args = parser.parse_args()

    if args.rt_barcode_file == "default":
        if args.level == "3":
            rtfile = RT3_FILE
            fix = 3
        else:
            rtfile = RT_FILE
            fix = 2
    else:
        rtfile = args.rt_barcode_file
        fix = 0
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
        if not line[0] in rtdict:
            sys.stderr.write("Sample sheet error at line " + str(line_num) + ". RT Barcode '" + line[0] + "' not valid.\n")
            error_flag = 1
        if not line[2] in genomes:
            sys.stderr.write("Sample sheet error at line " + str(line_num) + ". Reference Genome '" + line[2] + "' not valid.\n")
            error_flag = 1
        return error_flag

    def fix_line(line, fix):
        line = '%s\n' % ( line.strip() )
        line = line.split(",")
        line[1] = well_dict[line[0]]
        if fix == 0:
            return ",".join(line)
        if fix == 2:
            p1 = line[0].split("-")[0]
            p2 = line[0].split("-")[1]
            p1 = p1[0:2] + "{0:0=2d}".format(int(p1[2:]))
            p2 = p2[0:1] + "{0:0=2d}".format(int(p2[1:]))
            line[0] = p1 + "-" + p2
            return ",".join(line)
        if fix == 3:
            p1 = line[0].split("-")[0]
            p2 = line[0].split("-")[1]
            p1 = p1[0:1] + "{0:0=2d}".format(int(p1[1:]))
            p2 = p2[0:1] + "{0:0=2d}".format(int(p2[1:]))
            line[0] = p1 + "-" + p2
            return ",".join(line)

    def fix_line_exp(line, fix):
        line = '%s\n' % ( line.strip() )
        line = line.split(",")
        line[1] = well_dict[(line[0],line[3])]
        if fix == 0:
            return ",".join(line)
        if fix == 2:
            p1 = line[0].split("-")[0]
            p2 = line[0].split("-")[1]
            p1 = p1[0:2] + "{0:0=2d}".format(int(p1[2:]))
            p2 = p2[0:1] + "{0:0=2d}".format(int(p2[1:]))
            line[0] = p1 + "-" + p2
            return ",".join(line)
        if fix == 3:
            p1 = line[0].split("-")[0]
            p2 = line[0].split("-")[1]
            p1 = p1[0:1] + "{0:0=2d}".format(int(p1[1:]))
            p2 = p2[0:1] + "{0:0=2d}".format(int(p2[1:]))
            line[0] = p1 + "-" + p2
            return ",".join(line)

    # precount samples per well
    sample_dict = dict()
    well_dict = dict()
    with_exp = False
    with open(args.sample_sheet, 'r') as f:
        for line in f:
            line = line.split(",")
            if len(line) > 3:
                if line[1] in sample_dict:
                    sample_dict[line[1]].append((line[0],line[3]))
                else:
                    sample_dict[line[1]] = [(line[0],line[3])]
                well_dict[(line[0],line[3])] = line[1]
                with_exp = True
            elif with_exp:
                sys.stderr.write("One or more entries does not include Experiment label (third column) while others do.\n")
                sys.exit(11)
            else:
                if line[1] in sample_dict:
                    sample_dict[line[1]].append(line[0])
                else:
                    sample_dict[line[1]] = [line[0]]
                well_dict[line[0]] = line[1]

    for samp in sample_dict.keys():
        if len(sample_dict[samp]) > int(args.max_wells_per_sample):
            div = len(sample_dict[samp])/math.ceil(len(sample_dict[samp])/int(args.max_wells_per_sample))
            group_count = 1
            curr_count = 1
            for well in sample_dict[samp]:
                if curr_count <= div:
                    curr_count += 1  
                else:
                    curr_count = 1
                    group_count += 1
                well_dict[well] += "_fq_part" + str(group_count)

    sheet = open(args.sample_sheet)

    if with_exp:
        fix_func = fix_line_exp
    else:
        fix_func = fix_line

    # check for header
    topline_orig = sheet.readline()
    topline = topline_orig.strip().split(",")
    line_num = 1
    num_cols = len(topline)
    if num_cols == 3:
        if topline[1] =='Sample ID' and topline[2] == 'Reference Genome':
            sample_out = 'RT Barcode,Sample ID,Reference Genome\n'
        else:
            sample_out = 'RT Barcode,Sample ID,Reference Genome\n'
            topline_orig = fix_func(topline_orig, fix)
            check_line(topline_orig, line_num)
            sample_out = sample_out + topline_orig
    else:
        if topline[1] =='Sample ID' and topline[2] == 'Reference Genome':
            sample_out = topline_orig
        else:
            sample_out = 'RT Barcode,Sample ID,Reference Genome'
            for i in range(3, len(topline)):
                sample_out = sample_out + ",Column" + str(i)
            topline_orig = fix_func(topline_orig, fix)
            check_line(topline_orig, line_num)
            sample_out = sample_out + "\n" +  topline_orig
    # Check RT and Genomes against possible
    error_count = 0

    line_num = 1
    for line in sheet:
        line = fix_func(line, fix)
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

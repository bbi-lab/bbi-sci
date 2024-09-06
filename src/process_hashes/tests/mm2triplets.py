#!/usr/bin/env python3

import sys
from scipy.sparse import csc_matrix
from io import StringIO
from scipy.io import mmread
from scipy.sparse import coo_matrix


hash_sheet = '/net/bbi/vol1/data/regression_tests/sciRNAseq/reference_runs/nf.RNA3-072-a.ubuntu.ref_1/timecourse_hash2.txt'

key = sys.argv[1]

#
# Read MM file.
#
mtx_name = key + '.hashumis.mtx'
matrix = mmread(mtx_name)
matrix_coo = coo_matrix(matrix)

#
# Read row and column name files.
#
file_name = key + '.hashumis_hashes.txt'
fp = open(file_name, 'r')
row_names = []
for row_name in fp:
  row_names.append(row_name.strip())

file_name = key + '.hashumis_cells.txt'
fp = open(file_name, 'r')
col_names = []
for col_name in fp:
  col_names.append(col_name.strip())

#
# Read hash sheet file.
#
file_name = hash_sheet
fp = open(file_name, 'r')
hash_dict = {}
for line in fp:
  toks = line.split('\t')
  hash_name = toks[0].strip()
  hash_name = hash_name.replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')
  hash_seq  = toks[1].strip()
  hash_dict[hash_name] = hash_seq

for i, (r,c,d) in enumerate(zip(matrix_coo.row, matrix_coo.col, matrix_coo.data)):
  print("%s|%s|%s|%d" % (row_names[r], hash_dict[row_names[r]], col_names[c], d))

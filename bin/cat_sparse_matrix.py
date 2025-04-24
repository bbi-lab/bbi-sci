#!/usr/bin/env python3

#
# Concatenate sparse matrix files by column; that is,
# stack up cells given that the features are the same
# and in the same order in each of the files of
# triplets. Give the names of the files to concatenate
# on the command line. The names must end with
# 'matrix.mtx'.
#
# %%MatrixMarket matrix coordinate integer general
# %
# 60676 766 118308
# 8283 1 1
# 8912 1 1
#
# Notes:
#   o  all input matrix files must have exactly the
#      same row names and in the same order
#   o  all column names must be distinct
#   o  all matrices must be in COO format and have
#      integer values
#


import sys
import hashlib
import re
import argparse
import shutil ### avo


#
# Calculate the md5 checksum of a file.
#
def calculate_md5(file_path):
  hasher = hashlib.md5()
  with open(file_path, 'rb') as f:
    for chunk in iter(lambda: f.read(4096), b''):
        hasher.update(chunk)
  return hasher.hexdigest()


#
# Check that the files of row names (features)
# all have the same md5 checksum. We need the
# rows to be the same for all input matrices.
#
def check_rownames(file_name_list):
  md5_list = []
  num_file = 0
  for file_names in file_name_list:
    print(f'file names {file_names}')
    num_file += 1
    file_name = file_names[1]
    md5_file = calculate_md5(file_name)
    md5_list.append(md5_file)
  md5_1 = md5_list[0]
  result = md5_list and all(md5_list[0] == elem for elem in md5_list)
  print(md5_list)
  if(result == False):
    print('Error: features files differ.')
    sys.exit(-1)
#  print('Checked %d features files.' % (num_file))


#
# Get the matrix dimensions of all matrix files.
# The dimensions follow immediately the header comments.
# We use the dimensions to form the dimensions line of
# the output file, and to shift the column coordinates
# of the triplets.
#
def gather_matrix_dimensions(file_name_list):
  matrix_dimension_list = []
  num_file = 0
  for file_names in file_name_list:
    num_file += 1
    file_name = file_names[0]
    num_line = 0
    with open(file_name, 'r') as fh:
      for line in fh:
        num_line += 1
        if(num_line == 1):
          if(not re.match(r'%%MatrixMarket matrix coordinate (integer|real)', line)):
            print('Error: file %s has unexpected header line.' % (file_name), file=sys.stderr)
            print('  header line: %s' % (line))
            sys.exit(-1)
        elif(re.match(r'%', line)):
          continue
        else:
          toks = line.split()
#          print('%d %d %d' % (int(toks[0]), int(toks[1]), int(toks[2])))
          matrix_dimension_list.append([int(toks[0]), int(toks[1]), float(toks[2])])
          break
  return(matrix_dimension_list)          


#
# Concatenate the input matrices by adding columns to
# the output matrix file. The column values of the
# input matrix triplets are increased by the sum of
# the number of columns of all preceding input
# matrices. Each triplet is written to the output file
# immediately after shifting the column value.
#
# The cell name files are simply concatenated to form
# the output cell name file.
#
# The first feature name file is copied to the output
# feature name file.
#
def matrix_concatenation(file_name_list, matrix_dimension_list, out_rootname):
  out_matrix_name = '%s.hashumis.mtx' % (out_rootname)
  try:
    ofh = open(out_matrix_name, 'w')
  except:
    print('Error: unable to open output file \'%s\'' % (out_matrix_name), sys.stderr)
    sys.exit(-1)

  sum_col = 0
  sum_elem = 0
  for matrix_dim in matrix_dimension_list:
    sum_col += matrix_dim[1]
    sum_elem += matrix_dim[2]

  print('%%MatrixMarket matrix coordinate integer general', file=ofh)
#  print('%', file=ofh)
  print('%d %d %d' % (matrix_dimension_list[0][0], sum_col, sum_elem), file=ofh)

  #
  # It makes sense to process the comments and counts line
  # outside the main loop, when there is time to work on it.
  #
  cum_barcode_count = 0
  for ifile, file_names in enumerate(file_name_list):
    in_matrix_name = file_names[0]
#    print('==== %s' % (in_matrix_name), file=ofh)
    with open(in_matrix_name, 'r') as ifh:
      num_data = 0
      for line in ifh:
        if(re.match(r'%', line)):
          continue
        num_data += 1
        if(num_data == 1):
          continue
        toks = line.split()
        print('%s %d %s' % (toks[0], int(toks[1]) + cum_barcode_count, toks[2]), file=ofh)
    cum_barcode_count += matrix_dimension_list[ifile][1]

  #
  # Copy cell names to a file.
  #
  out_barcodes_name = '%s.hashumis_cells.txt' % (out_rootname)
  try:
    ofh = open(out_barcodes_name, 'w')
  except:
    print('Error: unable to open output file \'%s\'' % (out_barcodes_name), sys.stderr)
    sys.exit(-1)

  for ifile, file_names in enumerate(file_name_list):
    in_barcodes_name = file_names[2]
    num_lines = 0
    with open(in_barcodes_name, 'r') as ifh:
      for line in ifh:
        num_lines += 1
        print('%s' % (line.strip()), file=ofh)
    if(num_lines != matrix_dimension_list[ifile][1]):
      print('Error: inconsistent cell name count in file \'%s\' (%d != %d)' % (in_barcodes_name, num_lines, matrix_dimension_list[ifile][1]), file=sys.stderr)
      sys.exit(-1)

  ofh.close()

  #
  # Check for duplicate cell names.
  #
  with open(out_barcodes_name, 'r') as ifh:
    cell_name_set = set()
    for line in ifh:
      cell_name = line.strip()
      if(cell_name in cell_name_set):
        print('Error: duplicate cell names in count matrix.', file=sys.stderr)
        sys.exit(-1)
      cell_name_set.add(cell_name)

  #
  # Copy features to a file.
  #
  out_features_name = '%s.hashumis_hashes.txt' % (out_rootname)
  try:
    ofh = open(out_features_name, 'w')
  except:
    print('Error: unable to open output file \'%s\'' % (out_features_name), sys.stderr)
    sys.exit(-1)

  num_feature = 0
  in_features_name = file_names[1]
  with open(in_features_name, 'r') as ifh:
    for line in ifh:
      print('%s' % (line.strip()), file=ofh)
      num_feature += 1

  #
  # Check that the number of features is the same
  # as the number of matrix rows.
  #
  if(num_feature != matrix_dimension_list[0][0]):
    print('Error: number of row names != number of matrix rows.')
    sys.exit(-1)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to concatenate sparse matrix files, in triplet format, by column.')
  parser.add_argument('-i', '--input', required=True, default=None, nargs='+', help='Input sparse matrix filenames (required strings).')
  parser.add_argument('-o', '--output_root', required=False, default=None, help='Output files root name (required string).')
  parser.add_argument('-m', '--matrix_root', required=True, default=None, help='Matrix file root name (required string).')
  parser.add_argument('-f', '--feature_root', required=True, default=None, help='Features root name (required string).')
  parser.add_argument('-c', '--cell_root', required=True, default=None, help='Cells root name (required string).')
  parser.add_argument('-v', '--version', required=False, default=None, help='Write version string to stdout.')
  args = parser.parse_args()


  ### atv additions to account for only one hash file
  if len(args.input) == 1:
      print("Only one input file — forcing new writes to Nextflow outputs.")

      matrix_file = args.input[0]
      matrix_file_root = args.matrix_root
      feature_file_root = args.feature_root
      cell_file_root = args.cell_root
      path_base = matrix_file.replace(matrix_file_root, '')

      features_file = f'{path_base}{feature_file_root}'
      cells_file = f'{path_base}{cell_file_root}'

      out_rootname = args.output_root

      shutil.copyfile(matrix_file, f'{out_rootname}.hashumis.mtx')
      shutil.copyfile(features_file, f'{out_rootname}.hashumis_hashes.txt')
      shutil.copyfile(cells_file, f'{out_rootname}.hashumis_cells.txt')

      sys.exit(0)
  ### atv end


  file_name_list = []

  for i in range(len(args.input)):
    matrix_file_name = args.input[i]

    matrix_file_root = args.matrix_root
    features_file_root = args.feature_root
    cell_file_root = args.cell_root
    path_name = matrix_file_name.replace(matrix_file_root, '')

    matrix_name = matrix_file_name
    features_name = '%s%s' % (path_name, features_file_root)
    cells_name = '%s%s' % (path_name, cell_file_root)

    file_name_list.append([matrix_name, features_name, cells_name])
    print('matrix_name: %s' % (matrix_name))

  check_rownames(file_name_list)
  matrix_dimension_list = gather_matrix_dimensions(file_name_list)

  #
  # The output files have the names
  #   <out_rootname>.(raw|filtered).matrix.mtx
  #   <out_rootname>.(raw|filtered).features.tsv
  #   <out_rootname>.(raw|filtered).cells.tsv
  #
  out_rootname = args.output_root
  matrix_concatenation(file_name_list, matrix_dimension_list, out_rootname)


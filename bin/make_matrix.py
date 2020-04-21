#!/usr/bin/env python
import argparse
import sys
import numpy as np
from scipy.sparse import csc_matrix
import gzip

def read_old_sparse(filename, feature_dict):
    cell_dict = dict()
    row = list()
    col_names = list()
    col = list()
    data = list()
    col_count = 0
#    with open(filename, 'rt') as input_file:
    for line in filename:
        cell, gene, value = line.strip().split("\t")
        gene_index = feature_dict[gene]
        if cell in cell_dict:
            cell_index = cell_dict[cell]
        else:
            cell_index = col_count
            col_names.append(cell)
            cell_dict[cell] = col_count
            col_count += 1
        if col_count % 100000 == 0:
            print(col_count)
        row.append(gene_index)
        col.append(cell_index)
        data.append(int(value))
    return col,row,data,col_names
    
def write_list(l, filename):
    with open(filename, 'w') as output_file:
        for i in l:
            output_file.write('%s\n' % i)


def write_mtx_file(mat, column_names, key):
    from scipy.io import mmwrite

    write_list(list(column_names), key + ".cell_annotations.txt")
    
    mmwrite(target=key + ".umi_counts.mtx", a=mat)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to take umi rollup and generate matrix market sparse matrix files.')
    parser.add_argument('umi_rollup', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Text piped from std in for UMI rollup output file (cell, gene, value file.')
    parser.add_argument('--gene_annotation', required=True, help='Gene annotation file with all genes to be included as first column.')
    parser.add_argument('--key', help='sample name for output file.')
    args = parser.parse_args()

    feature_dict = dict()
    count = 0
    with open(args.gene_annotation, 'r') as gene_file:
        for line in gene_file:
            gene_id = line.strip().split("\t")[0]
            feature_dict[gene_id] = count
            count += 1
    col,row,data,col_names = read_old_sparse(args.umi_rollup, feature_dict)
    
    sparse_mat = csc_matrix((np.array(data), (np.array(row), np.array(col))), shape=(len(feature_dict), len(col_names)), dtype='int32')
    write_mtx_file(sparse_mat, col_names, args.key)



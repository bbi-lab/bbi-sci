#!/usr/bin/env python
import barcodeutils as bu
import argparse
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
from scipy.sparse import csc_matrix

def load_hash_layout(file_path):
    """
    Function that loads the sample layout file to a hash lookup table.
    """
    lookup = {}
    fh = open(file_path)
    for line in fh:
        entry = line.strip().split("\t")
        lookup[entry[1]] = entry[0].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')
    fh.close()
    return lookup

def write_list(l, filename):
    with open(filename, 'w') as output_file:
        for i in l:
            output_file.write('%s\n' % i)


def write_mtx_file(mat, column_names, row_names, key):
    from scipy.io import mmwrite

    write_list(list(column_names), key + ".hashumis_cells.txt")
    write_list(list(row_names), key + ".hashumis_hashes.txt")

    mmwrite(target=key + ".hashumis.mtx", a=mat)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to combine R1 and R2 file from sci run into a single file with barcodes in the fastq header.')
    parser.add_argument('--hash_sheet', required=True, help='Path to hash sample sheet.')
    parser.add_argument('--fastq', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for fastq.')
    parser.add_argument('--hash_edit_distance', type=int, required=False, default=1, help='Allowed edit distance for hashes.')
    parser.add_argument('--key', required=True, help='Key for file name for output file.')
    args = parser.parse_args()

    hash_lookup = load_hash_layout(args.hash_sheet)
    hash_whitelist = bu.construct_mismatch_to_whitelist_map(hash_lookup, edit_distance = args.hash_edit_distance)

    # Dictionary to count hash reads per cell
    hashdict = dict()
    # Dictionary to count hash UMIs per cell
    hashumi = dict()
    for h in hash_lookup.keys():
        hashdict[h] = dict()
        hashumi[h] = dict()

    cells = set()

    # Dictionary to count total hashes assigned
    hashcounts = dict()
    for i in range(args.hash_edit_distance + 1):
        hashcounts[i] = 0
    
    inp_handle = FastqGeneralIterator(args.fastq)
    for r in inp_handle:
        line1, line2, line4 = r
        is_hash = False
        hashbc = line2[0:10]
        polya = line2[11:15]
        for i in range(args.hash_edit_distance + 1):
            if hashbc in hash_whitelist[i] and polya == line2[11:15] == "AAAA":
                hashval = hash_whitelist[i][hashbc]
                is_hash = True
                break
        if not is_hash:
            continue
        cell_barc = line1.split("|")[2:-1]
        cell_barc = "_".join(cell_barc)
        umi = line1.split("|")[-1]
        cells.add(cell_barc)
        if cell_barc in hashdict[hashval]:
            lpre = len(hashdict[hashval][cell_barc])
            hashdict[hashval][cell_barc].add(umi)
            lpost = len(hashdict[hashval][cell_barc])
            if lpost > lpre:
                hashcounts[i] += 1
        else:
            hashdict[hashval][cell_barc] = set(umi)
            hashcounts[i] += 1
                            
    data = []
    row = []
    col = []
    hashname = []
    for c,cell in enumerate(cells): 
        for h,hashkey in enumerate(hashdict.keys()):
            if cell in hashdict[hashkey]:
                col.append(c)
                row.append(h)
                hashname.append(hash_lookup[hashkey])
                data.append(len(hashdict[hashkey][cell]))

    hashname = []
    for hashkey in hashdict.keys():
        hashname.append(hash_lookup[hashkey])

    sparse_mat = csc_matrix((np.array(data), (np.array(row), np.array(col))), shape=(len(hashdict), len(cells)), dtype='int32')
    write_mtx_file(sparse_mat, cells, hashname, args.key)

    log = open(args.key + "_hash.log", 'w')
    for i in range(args.hash_edit_distance + 1):
        log.write("Hash UMIs detected with " + str(i) + " correction: " + str(hashcounts[i]) + "\n")
    log.close()

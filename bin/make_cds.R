#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparse)
    library(Matrix)
    library(monocle3)
})

parser = argparse::ArgumentParser(description='Script to make final cds per sample.')
parser$add_argument('matrix', help='File of umi count matrix.')
parser$add_argument('cell_data', help='File of cell data.')
parser$add_argument('gene_data', help='File of gene data.')
args = parser$parse_args()

sample_name <- stringr::str_split_fixed(args$matrix, "\\.txt\\.", 2)[,1]

cds <- load_mtx_data(mat_path = args$matrix, gene_anno_path = args$gene_data, cell_anno_path = args$cell_data, umi_cutoff=100)

saveRDS(cds, file=paste0(sample_name, "_cds.RDS"))

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparse)
    library(Matrix)
    library(monocle3)
    library(DropletUtils)
})

parser = argparse::ArgumentParser(description='Script to make final cds per sample.')
parser$add_argument('matrix', help='File of umi count matrix.')
parser$add_argument('cell_data', help='File of cell data.')
parser$add_argument('gene_data', help='File of gene data.')
parser$add_argument('key', help='The sample name prefix.')
parser$add_argument('output_file', help='The output filename.')
args = parser$parse_args()

sample_name <- args$key

ed_umi_cutoff <- 0
cds <- load_mm_data(mat_path = args$matrix, feature_anno_path = args$gene_data, 
                    cell_anno_path = args$cell_data, umi_cutoff=ed_umi_cutoff,
                    feature_metadata_column_names=c('gene_short_name'), sep="")

#
# Notes:
#   o monocle3::load_mm_data() includes cells with >= umi_cutoff=100
#   o emptyDrops excludes cells with <= 'ignore' UMI
#
ed_lower <-        60
ed_niters <-       10000
ed_test_ambient <- TRUE
ed_ignore <-       99
ed_alpha <-        NULL
ed_round <-        TRUE
ed_retain <-       Inf

emptyDrops_out <- emptyDrops(m=SingleCellExperiment::counts(cds),
                             lower=ed_lower,
                             niters=ed_niters,
                             test.ambient=ed_test_ambient,
                             ignore=ed_ignore, 
                             alpha=ed_alpha,
                             round=ed_round,
                             retain=ed_retain)

metadata(emptyDrops_out)$ignore <-  ed_ignore
metadata(emptyDrops_out)$round <-   ed_round

saveRDS(object=emptyDrops_out, file=args$output_file)



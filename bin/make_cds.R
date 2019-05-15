#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparse)
    library(Matrix)
    library(monocle)
})

parser = argparse::ArgumentParser(description='Script to make final cds per sample.')
parser$add_argument('matrix', help='File of umi count matrix.')
parser$add_argument('cell_data', help='File of cell data.')
parser$add_argument('gene_data', help='File of gene data.')
args = parser$parse_args()

get_mat <- function(mat.path, gene.annotation.path, cell.annotation.path) {
    df = read.table(
        mat.path,
        col.names = c("gene.idx", "cell.idx", "count"),
        colClasses = c("integer", "integer", "integer"))

    gene.annotations = read.table(
        gene.annotation.path,
        col.names = c("id", "gene_short_name"),
        colClasses = c("character", "character"))

    cell.annotations = read.table(
        cell.annotation.path,
        col.names = c("cell"),
        colClasses = c("character"))

    rownames(gene.annotations) = gene.annotations$id
    rownames(cell.annotations) = cell.annotations$cell

    df = rbind(df, data.frame(
        gene.idx = c(1, nrow(gene.annotations)),
        cell.idx = rep(nrow(cell.annotations)+1, 2),
        count = c(1, 1)))

    mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
    mat = mat[, 0:(ncol(mat)-1), drop=FALSE]

    rownames(mat) = gene.annotations$id
    colnames(mat) = cell.annotations$cell
    mat
}

mat <- get_mat(args$matrix, args$gene_data, args$cell_data)

cell.annotations <- data.frame(cell = colnames(mat))
row.names(cell.annotations) <- cell.annotations$cell

gene.annotations <- data.frame(gene = row.names(mat))
row.names(gene.annotations) <- gene.annotations$gene

pd = new("AnnotatedDataFrame", data = cell.annotations)
fd = new("AnnotatedDataFrame", data = gene.annotations)
cds = newCellDataSet(mat, phenoData = pd, featureData = fd, expressionFamily = VGAM::negbinomial.size())
pData(cds)$n.umi = Matrix::colSums(exprs(cds))

sample_name <- stringr::str_split_fixed(args$matrix, "\\.\\.", 2)[,1]

saveRDS(cds, file=paste0(sample_name, "_cds.RDS"))

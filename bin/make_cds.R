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
parser$add_argument('gene_bed', help='Bed file of gene info.')
parser$add_argument('key', help='The sample name prefix.')
parser$add_argument('umi_cutoff', help='UMI cutoff to count as a cell.')
args = parser$parse_args()

sample_name <- args$key

cds <- load_mm_data(mat_path = args$matrix, feature_anno_path = args$gene_data, 
                    cell_anno_path = args$cell_data, umi_cutoff=args.umi_cutoff,
                    feature_metadata_column_names=c('gene_short_name'), sep="")
gene_bed <- read.table(args$gene_bed)
row.names(gene_bed) <- gene_bed$V4
names(gene_bed) <- c("chromosome", "bp1", "bp2", "id", "x", "strand")

temp <- gene_bed[row.names(fData(cds)),]
fData(cds)[,c("id", "chromosome", "bp1", "bp2", "gene_strand")] <- temp[,c("id", "chromosome", "bp1", "bp2", "strand")]
#fData(cds)$gene_biotype <- gene_bed[row.names(fData(cds)),"gene_biotype"]

mt <- row.names(fData(cds)[!is.na(fData(cds)$chromosome) & (fData(cds)$chromosome %in% c("MT", "MtDNA", "Mt", "HUMAN_MT", "MOUSE_MT", "mitochondrion_genome")),])
mt_cds <- cds[mt,]
pData(cds)$perc_mitochondrial_umis <- Matrix::colSums(exprs(mt_cds))/Matrix::colSums(exprs(cds)) * 100

#rt <- row.names(fData(cds)[!is.na(fData(cds)$gene_biotype) & fData(cds)$gene_biotype == "rRNA",])
#rt_cds <- cds[rt,]
#pData(cds)$perc_rRNA_umis <- Matrix::colSums(exprs(rt_cds))/Matrix::colSums(exprs(cds)) * 100

qc <- as.data.frame(pData(cds))[,c("cell", "n.umi", "perc_mitochondrial_umis")]
write.csv(qc, file=paste0(sample_name, "_cell_qc.csv"), quote=FALSE, row.names = FALSE)

writeMM(exprs(cds), paste0(sample_name, "_for_scrub.mtx"))

saveRDS(cds, file=paste0(sample_name, "_cds.RDS"))

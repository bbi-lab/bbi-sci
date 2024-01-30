#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(monocle3)
  library(data.table)
  library(pbapply)
})

parser = argparse::ArgumentParser(description='Script to assign hash per sample')
parser$add_argument('key', help='Sample name')
parser$add_argument('hash_matrix', help='File of hash umi count matrix')
parser$add_argument('cell_list', help='File with list of cell names with hash umis')
parser$add_argument('hash_list', help='File with list of hash names')
parser$add_argument('cds', help='cds object in RDS format')
parser$add_argument('umis_per_cell', help='File with list of umis per cell barcode')
parser$add_argument('hash_umi_cutoff', help='min number of hash umis to filter out')
args = parser$parse_args()

# Takes in cell hash matrix and background hash frequencies (determined by hash umi cutoff)
# and find pvalues of hash frequencies
chisq_vs_background <- function(test_hash_matrix, hash_frequencies){
  hash_frequencies_nz = which(hash_frequencies > 0)
  hash_frequencies = hash_frequencies[hash_frequencies_nz]
  pvals= pbapply(test_hash_matrix[,hash_frequencies_nz], 1, function(x) {
    tryCatch({
      res = chisq.test(x, p=hash_frequencies,  simulate.p.value = FALSE)
      unlist(res[["p.value"]])
    }, error = function(e) { 1.0 }
    )
  })
  return(pvals)
}

# Takes in a hash sparse matrix, cells from cds object and background cell hashes
# (determined by cells with less than hash umi cutoff) to assign hash labels 
# and determine the top to second best hash oligo.
# Creates a corrected hash table 
assign_hash_labels <- 
  function(hash_matrix,
           test_cell_hashes, 
           background_cell_hashes, 
           qval_thresh = 0.05, 
           downsample_rate=NULL){
    background_hash_matrix = hash_matrix[rownames(hash_matrix) %in% background_cell_hashes,]
    background_hash_matrix = background_hash_matrix[rowSums(background_hash_matrix)>0,]

    background_hash_frequencies = colSums(background_hash_matrix)/sum(colSums(background_hash_matrix))
    
    
    test_hash_matrix = hash_matrix[rownames(hash_matrix) %in% test_cell_hashes,]
    
    if (is.null(downsample_rate) == FALSE){
      test_hash_matrix = floor(test_hash_matrix * downsample_rate)
    }
    pvals = 
      chisq_vs_background(test_hash_matrix, 
                          hash_frequencies = background_hash_frequencies)
    qvals = p.adjust(pvals)

    expected_background_hashes = outer(rowSums(test_hash_matrix), background_hash_frequencies)

    background_subtracted_test_hashes = test_hash_matrix - expected_background_hashes

    background_subtracted_test_hashes[background_subtracted_test_hashes < 0] = 0


    #hash_hits = background_subtracted_test_hashes * (qvals < qval_thresh)
    top_to_second_best_ratios = apply(background_subtracted_test_hashes, 1, function(x) { y = sort(x, decreasing=TRUE); y[1] / y[2]})
    
    top_hash = apply(background_subtracted_test_hashes, 1, function(x) {
      m = which(x == max(x));
      if (length(m) > 1)
        return (NA)
      colnames(background_subtracted_test_hashes)[which(x == max(x))]
    })
    #unambiguous_hits = which(top_to_second_best_ratios > min_best_vs_second_best_ratio)

    hash_label_df = data.frame(cell = rownames(test_hash_matrix),
                               hash_umis = rowSums(test_hash_matrix),
                               pval = pvals,
                               qval = qvals,
                               #hit = qvals < qval_thresh,
                               top_to_second_best_ratio = top_to_second_best_ratios,
                               #unambiguous_hit = top_to_second_best_ratios > min_best_vs_second_best_ratio,
                               top_oligo = top_hash)
  }

cds <- readRDS(args$cds)

# Extract meta info from cell name 
df <- as.data.frame(colData(cds))
meta_types <-  meta_types <- c("P5_barcode", "P7_barcode", "RT_barcode", "Ligation_barcode")
meta<- separate(df, cell, into=meta_types, sep="_", remove=FALSE)

# Fill in barcode information in cds object from meta object 
for (m in meta_types) {
  colData(cds)[,m] <- meta[[m]]
}

# Extract pcr plate number
cds$plate <- sapply(strsplit(as.character(meta$RT_barcode), "-"), `[`, 1)

# Get list of cell names with hash umis 

cell_list <- fread(args$cell_list,
                  header = FALSE,
                  data.table = F)[,1]

hash_list = fread(args$hash_list, 
                  header = FALSE, 
                  data.table = F)[,1]


# Hash sparse matrix with hash umi counts for each cell 
# Add hash names to rows and cell names with hashes to column 
# Effectively, a hash by cell matrix 

hash_mtx = Matrix::readMM(args$hash_matrix)
hash_mtx = as(hash_mtx, "dgCMatrix")
rownames(hash_mtx) = hash_list
colnames(hash_mtx) = cell_list

# Transpose hash matrix so hashes are columns and cells are rows
hash_mtx = t(hash_mtx)

# Grab sci umis 
rna_umis = fread(args$umis_per_cell,
                 header = F, data.table = F, col.names = c("cell", "n.umi"))

# Filter for umis that are less than 5
# Determine background cell hashes to find top_to_second_best_ratio
background_cell_hashes =
  rna_umis$cell[rna_umis$n.umi < args$hash_umi_cutoff] %>%
  as.character()

# Assign number of hash umis to each cell, pvals, qvals, top_to_second_best_ration, and the top oligo
corrected_hash_table = assign_hash_labels(hash_matrix = hash_mtx,
                                          test_cell_hashes = cds$cell,
                                          background_cell_hashes = background_cell_hashes)

sample_name = args$key
fwrite(corrected_hash_table, file=paste0(sample_name, "_hash_table.csv"), sep = ",")

# merge hash table with cds to assign to cells
merged = as.data.table(merge(x=corrected_hash_table, y=colData(cds), by = "cell",all.x=FALSE, all.y=TRUE))
merged <- merged %>% mutate_at(vars("hash_umis"), ~replace_na(., 0)) # add 0 if a cell has no hash

# suppressMessages(merged <- splitstackshape::cSplit(merged, "top_oligo", "."))

colData(cds)$hash_umis = merged$hash_umis
colData(cds)$pval = merged$pval
colData(cds)$qval = merged$qval
colData(cds)$top_to_second_best_ratio = merged$top_to_second_best_ratio
colData(cds)$top_oligo = merged$top_oligo

# Drop any cells with less than hash umi cutof
cds_filt <- cds[,colData(cds)$hash_umis >= args$hash_umi_cutoff]
# cds_filt <- cds_filt[,colData(cds_filt)$top_to_second_best_ratio > args$hash_ratio]

saveRDS(cds_filt,file=paste0(sample_name, "_cds.RDS"))

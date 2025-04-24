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
parser$add_argument('hash_umi_cutoff', type='integer', help='min number of hash umis to determine top to second best ratio')
parser$add_argument('hash_ratio', help='min top to second best hash ratio. Default is false and not filtered')
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
    }, error = function(e) { 1.0
     }
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
    best_hash_umi = apply(background_subtracted_test_hashes, 1, function(x) { y = sort(x, decreasing=TRUE); round(y[1],4)})
    second_best_hash_umi = apply(background_subtracted_test_hashes, 1, function(x) { y = sort(x, decreasing=TRUE); round(y[2],4)})
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
                               top_oligo = top_hash,
                               best_hash_umi = best_hash_umi,
                               second_best_hash_umi = second_best_hash_umi)
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

# Extract RT plate number
cds$RT_plate <- sapply(strsplit(as.character(meta$RT_barcode), "-"), `[`, 1)

if ((file.info(args$cell_list)$size==0)) {
  saveRDS(cds,file=paste0(args$key, "_cds.RDS"))
  file.create(paste0(args$key, "_hash_table.csv"))
  quit(save = "no", status = 0)
}

# Read txt file of cell names with hash umis 
cell_list <- fread(args$cell_list,
                  header = FALSE,
                  data.table = F)[,1]

# Read txt file with hash oligo names 
hash_list = fread(args$hash_list, 
                  header = FALSE, 
                  data.table = F)[,1]

# Hash sparse matrix with hash umi counts for each cell 
# Add hash names to rows and cell names with hashes to column 
# Effectively, a hash by cell matrix 

hash_mtx = Matrix::readMM(args$hash_matrix)
hash_mtx = as(hash_mtx, "dgCMatrix")

# Set hash oligo names as rows and cell names as columns
rownames(hash_mtx) = hash_list
colnames(hash_mtx) = cell_list

# Transpose hash matrix so hashes are columns and cells are rows
hash_mtx = t(hash_mtx)

# Grab sci umis 
rna_umis = fread(args$umis_per_cell,
                 header = F, data.table = F, col.names = c("cell", "n.umi"))

# Filter for umis that are less than specified cutoff 
# Determine background cell hashes to find top_to_second_best_ratio
background_cell_hashes =
  rna_umis$cell[rna_umis$n.umi < args$hash_umi_cutoff] %>%
  as.character()

# Assign number of hash umis to each cell, pvals, qvals, top_to_second_best_ratio, and the top oligo
corrected_hash_table = assign_hash_labels(hash_matrix = hash_mtx,
                                          test_cell_hashes = cds$cell,
                                          background_cell_hashes = background_cell_hashes)

sample_name = args$key

# Drop any cells with less than hash umi cutoff
# corrected_hash_table <- filter(corrected_hash_table, hash_umis >= args$hash_umi_cutoff)
fwrite(corrected_hash_table, file=paste0(sample_name, "_hash_table.csv"), sep = ",")

# merge hash table with cds to assign to cells
# If hash table is empty, put NA for each hash column 
if (dim(corrected_hash_table)[1] != 0) {
  # merged = as.data.table(merge(x=corrected_hash_table, y=colData(cds), by = "cell",all.x=FALSE, all.y=TRUE))
  # merged <- merged %>% mutate_at(vars("hash_umis"), ~replace_na(., 0)) # add 0 if a cell has no hash

  # suppressMessages(merged <- splitstackshape::cSplit(merged, "top_oligo", "."))

  # hash_df <- merged %>% 
  #   select("cell", "hash_umis", "pval", "qval", "top_to_second_best_ratio", 
  #          "top_oligo", "best_hash_umi", "second_best_hash_umi")
  # 
  # to_merge <- data.frame(colData(cds)) %>% left_join(hash_df, by = "cell")
  # colData(cds) <- as(to_merge, "DataFrame")

  
  corrected_hash_table <- corrected_hash_table[,c('cell', 'hash_umis', 'pval', 'qval',
                                                  'top_to_second_best_ratio', 'top_oligo', 
                                                  'best_hash_umi', 'second_best_hash_umi')]
  
  # Bind the hash table columns to the colData(cds) keeping the colData(cds) dimension and row order.
  col_data_merged = as.data.frame(left_join(x=as.data.frame(colData(cds)), y=corrected_hash_table, by = "cell", keep=FALSE))
  # Add 0 if a cell has no hash
  col_data_merged <- col_data_merged %>% mutate_at(vars("hash_umis"), ~replace_na(., 0)) 
  colData(cds) <- as(col_data_merged, "DataFrame")

  # Drop any cells with less than top to second best hash ratio cutoff if args$hash_ratio is not false 
  if ("false" != "false") {
    cds  <- cds[,!is.na(colData(cds)$top_to_second_best_ratio) & colData(cds)$top_to_second_best_ratio >= 5 ]
  } 
}

saveRDS(cds,file=paste0(sample_name, "_cds.RDS"))

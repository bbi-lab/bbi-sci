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
parser$add_argument('--upper_umi_cutoff', type='integer', help='Optional. Specifies an upper UMI cutoff. Default is no upper UMI cutoff')
parser$add_argument('--hash_ratio', type='double', help='Optional. Specifies ratio of hash UMIs to determine top hash oligo')
parser$add_argument('--hash_umi_cutoff', type='integer', help='Optional. Specifies hash umi cutoff')
args = parser$parse_args()

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


assign_hash_labels <- 
  function(hash_matrix,
           test_cell_hashes, 
           background_cell_hashes, 
           qval_thresh = 0.05, 
           downsample_rate=NULL){
    background_hash_matrix = hash_matrix[rownames(hash_matrix) %in% background_cell_hashes,]
    background_hash_matrix = background_hash_matrix[rowSums(background_hash_matrix)>0,]
    

    # print("background_cell_hashes")
    # print(head(background_cell_hashes))
    # print("background_hash_matrix")
    # print(head(background_hash_matrix))

    background_hash_frequencies = colSums(background_hash_matrix)/sum(colSums(background_hash_matrix))
    
    
    test_hash_matrix = hash_matrix[rownames(hash_matrix) %in% test_cell_hashes,]
    
    if (is.null(downsample_rate) == FALSE){
      test_hash_matrix = floor(test_hash_matrix * downsample_rate)
    }
    pvals = 
      chisq_vs_background(test_hash_matrix, 
                          hash_frequencies = background_hash_frequencies)
    qvals = p.adjust(pvals)
    
    # print("background_hash_frequencies")
    # print(background_hash_frequencies)
    # print("test_hash_matrix")
    # print(head(test_hash_matrix))

    expected_background_hashes = outer(rowSums(test_hash_matrix), background_hash_frequencies)

    # print("expected_background_hashes")
    # print(head(expected_background_hashes))

    background_subtracted_test_hashes = test_hash_matrix - expected_background_hashes

    # print("background sub test hashs 1")
    # print(head(background_subtracted_test_hashes))
    # print(dim(background_subtracted_test_hashes))

    background_subtracted_test_hashes[background_subtracted_test_hashes < 0] = 0



    # print("background sub test hashs 2")
    # print(head(background_subtracted_test_hashes))
    # print(dim(background_subtracted_test_hashes))


    #hash_hits = background_subtracted_test_hashes * (qvals < qval_thresh)
    top_to_second_best_ratios = apply(background_subtracted_test_hashes, 1, function(x) { y = sort(x, decreasing=TRUE); y[1] / y[2]})
    
    top_hash = apply(background_subtracted_test_hashes, 1, function(x) {
      m = which(x == max(x));
      if (length(m) > 1)
        return (NA)
      colnames(background_subtracted_test_hashes)[which(x == max(x))]
    })
    #unambiguous_hits = which(top_to_second_best_ratios > min_best_vs_second_best_ratio)

    # print("rownames of test hash matrix")
    # print(head(rownames(test_hash_matrix)))
    # print(length(rownames(test_hash_matrix)))
    # print("rowSumsof test hash matrix")
    # print(head(rowSums(test_hash_matrix)))
    # print(length(rowSums(test_hash_matrix)))

    # print("pvals")
    # print(head(pvals))
    # print(length(pvals))
    # print("qvals")
    # print(head(qvals))
    # print(length(qvals))
    print("top to second")
    print(head(top_to_second_best_ratios))
    print(dim(top_to_second_best_ratios))
    # print("top hash")
    # print(head(top_hash))
    # print(length(top_hash))

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

# summary(cds$n.umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1       1       2      60       4 1040674 

##### This step splits barcode information out into:
##### pcr_plate, rt_plate_well, lig_well
##### could just use the "extract_barcode" function in our pipeline
##### Doesn't seem this information is used currently. Can just leave it for the wrapped data
cell_meta_data = colData(cds) %>%
  as_tibble() %>%
  select(cell) %>%
  separate(cell, into = c("pcr_p5", "pcr_p7", "rt_plate_well", "lig_well"), sep = "_") %>%
  mutate(pcr_plate = paste(str_sub(pcr_p7, start = 1, end = 1), str_sub(pcr_p5, start = 2, end = 3), sep = "")) %>%
  select(pcr_plate, everything(), -pcr_p5, -pcr_p7)
colData(cds)$pcr_plate = cell_meta_data$pcr_plate
colData(cds)$rt_plate_well = cell_meta_data$rt_plate_well
colData(cds)$lig_well = cell_meta_data$lig_well

# DataFrame with 6 rows and 7 columns
#                                     cell Size_Factor     n.umi perc_mitochondrial_umis   pcr_plate rt_plate_well    lig_well
# <character>   <numeric> <numeric>               <numeric> <character>   <character> <character>
#   A05_D01_P01-A01_LIG1     A05_D01_P01-A01_LIG1    0.952466         3                       0         D05       P01-A01        LIG1
# A05_D01_P01-A01_LIG10   A05_D01_P01-A01_LIG10    0.634977         2                       0         D05       P01-A01       LIG10
# A05_D01_P01-A01_LIG100 A05_D01_P01-A01_LIG100    0.317489         1                       0         D05       P01-A01      LIG100
# A05_D01_P01-A01_LIG102 A05_D01_P01-A01_LIG102    0.634977         2                       0         D05       P01-A01      LIG102
# A05_D01_P01-A01_LIG103 A05_D01_P01-A01_LIG103    1.587443         5                       0         D05       P01-A01      LIG103
# A05_D01_P01-A01_LIG104 A05_D01_P01-A01_LIG104    0.317489         1                       0         D05       P01-A01      LIG104

# UMI cutoff
# Prob don't need the lower cutoff if running through our pipeline 
lower_cutoff = 100 #### need to think about if umi cutoff for cds is not 100
upper_cutoff = args$upper_umi_cutoff

#### should already have a lower cutoff of 100 so may not need that portion 
cds = cds[,colData(cds)$n.umi >= lower_cutoff & colData(cds)$n.umi <= upper_cutoff]
# ncol(cds) # number cells left after filtering 


# Gets list of cell names with hash umis 
# cell_list = fread("CHEMFISH.ATOH7/CHEMFISH.ATOH7.hashumis_cells.txt", 
#                   header = FALSE, 
#                   data.table = F)[,1]

cell_list <- fread(args$cell_list,
                  header = FALSE,
                  data.table = F)[,1]

# head(cell_list)
# "A05_D04_P01-E11_LIG19"  "E05_D10_P01-D07_LIG99"  "A06_E04_P01-A06_LIG19"  "A06_E10_P01-C12_LIG42"  "C06_E04_P01-H12_LIG79"  "D05_D05_P01-B02_LIG131"

# length(cell_list)
# 1397012


hash_list = fread(args$hash_list, 
                  header = FALSE, 
                  data.table = F)[,1]
                
# head(hash_list)
# length(hash_list)

# WT.54.hpf.P02.A1" "WT.54.hpf.P02.A2" "WT.54.hpf.P02.A3" "WT.54.hpf.P02.A4" "WT.54.hpf.P02.A5" "WT.54.hpf.P02.A6" 
# 148 

# Hash sparse matrix with hash umi counts for each cell 
# Add hash names to rows and cell names with hashes to column 
# Effectively, a hash by cell matrix 

hash_mtx = Matrix::readMM(args$hash_matrix)
hash_mtx = as(hash_mtx, "dgCMatrix")
# hash_mtx[1:5, 1:5]
# dim(hash_mtx)
rownames(hash_mtx) = hash_list
colnames(hash_mtx) = cell_list

print("hash_mtx")
print(head(hash_mtx))

# Flip it so that the hashes are columns and the cells are rows
hash_mtx = t(hash_mtx)
# Comes in as a sparse matrix
# hash_mtx[1:5, 1:5]

# Grab sci umis 
#### This step seems redundant as this information is in the cds 
# rna_umis = fread("CHEMFISH.ATOH7/umis_per_cell_barcode.txt", 
#                  header = F, data.table = F, col.names = c("Cell", "n.umi"))

# rna_umis_test <- data.frame(Cell=colData(cds)$cell, n.umi=colData(cds)$n.umi)
# dim(rna_umis)
# # dim(rna_umis)
# # [1] 3209903       2
# 
# 
# rna_umis %>% head
# Cell n.umi
# 1   A05_D01_P01-A01_LIG1     3
# 2  A05_D01_P01-A01_LIG10     2
# 3 A05_D01_P01-A01_LIG100     1
# 4 A05_D01_P01-A01_LIG102     2
# 5 A05_D01_P01-A01_LIG103     5
# 6 A05_D01_P01-A01_LIG104     1



rna_umis = fread(args$umis_per_cell,
                 header = F, data.table = F, col.names = c("cell", "n.umi"))

# rna_umis_test <- data.frame(cell=colData(cds)$cell, n.umi=colData(cds)$n.umi)

# Filter for umis that are less than 5
#### Determine background cell hashes to find top_to_second_best_ratio
background_cell_hashes =
  rna_umis$cell[rna_umis$n.umi < 5] %>%
  as.character()



# Grabbing cell names 
# test_cells =
#   colnames(cds) %>%
#   as.character()


#### Assigns number of hash umis to each cell, pvals, qvals, top_to_second_best_ration, and the top oligo

# print("hash mtx")
# print(dim(hash_mtx))
# print("cds$cell")
# print(head(cds$cell))
# print(length(cds$cell))
corrected_hash_table = assign_hash_labels(hash_matrix = hash_mtx,
                                          test_cell_hashes = cds$cell,
                                          background_cell_hashes = background_cell_hashes)


# cell hash_umis          pval          qval top_to_second_best_ratio              top_oligo
# A05_D04_P01-E11_LIG19 A05_D04_P01-E11_LIG19        13 2.135124e-143 2.684769e-138                 3.361410       SB.72.28.1.P1.H4
# A06_E04_P01-A06_LIG19 A06_E04_P01-A06_LIG19        30  4.420273e-87  3.572376e-82                 4.654625 ATOH7.KO.72.hpf.P02.D3
# A06_E10_P01-C12_LIG42 A06_E10_P01-C12_LIG42        19  1.037784e-14  8.152829e-11                 1.965502 ATOH7.KO.72.hpf.P02.D3
# C06_E04_P01-H12_LIG79 C06_E04_P01-H12_LIG79       319  0.000000e+00  0.000000e+00                 1.211253 ATOH7.KO.54.hpf.P02.B3
# F06_E06_P01-G12_LIG3   F06_E06_P01-G12_LIG3        95  0.000000e+00  0.000000e+00                21.076158 ATOH7.KO.72.hpf.P02.D1
# E05_D10_P01-C12_LIG63 E05_D10_P01-C12_LIG63        10 9.661131e-164 1.331149e-158                      Inf       LY.72.28.2.P1.D8


# Save this 
sample_name = args$key
print("sample_name")
print(sample_name)
fwrite(corrected_hash_table, file=paste0(sample_name, "_hash_table.csv"), sep = ",")

# merge hash table with cds to assign to cells
#### not sure if we need this step. seems to just add another cell name column that already exists 
#### Looks like the cds cell is labeled with "cell" and hash table is labeled as "Cell" 
#### Adjusted for this in the source script to be consistent with cds label and can remove this step
# corrected_hash_table$cell = corrected_hash_table$Cell 

merged = as.data.table(merge(x=corrected_hash_table, y=colData(cds), by = "cell",all.x=FALSE, all.y=TRUE))

#                         cell hash_umis   pval          qval           top_to_second_best_ratio      top_oligo Size_Factor n.umi perc_mitochondrial_umis P5_barcode
# 1: A05_D01_P01-A01_LIG107        11  5.218253e-29  1.002896e-24                 1.718059 ATOH7.KO.72.hpf.P02.D8    37.46366   118               0.8474576        A05
# 2: A05_D01_P01-A01_LIG160        38 6.462316e-289 1.195186e-283                11.537028       WT.72.hpf.P02.C3   504.80688  1590               5.4716981        A05
# 3: A05_D01_P01-A01_LIG163        20  0.000000e+00  0.000000e+00                 3.302363   WntC59.72.28.7.P1.A7   125.72549   396               3.2828283        A05
# 4:  A05_D01_P01-A01_LIG19         7  1.477237e-10  7.544251e-07                 2.496314 ATOH7.KO.72.hpf.P02.D4    59.05288   186               0.0000000        A05
# 5:  A05_D01_P01-A01_LIG23        11  2.701848e-82  2.056836e-77                 1.213414     DMSO.72.28.5.P1.G5   103.81877   327               0.0000000        A05
# 6:  A05_D01_P01-A01_LIG39        35 3.664147e-180 5.365227e-175                 1.741385   TGFB2.72.28.0.P1.F11   150.80709   475               0.4210526        A05
# P7_barcode RT_barcode Ligation_barcode plate num_genes_expressed
# 1:        D01    P01-A01           LIG107   P01                 101
# 2:        D01    P01-A01           LIG160   P01                 878
# 3:        D01    P01-A01           LIG163   P01                 274
# 4:        D01    P01-A01            LIG19   P01                 151
# 5:        D01    P01-A01            LIG23   P01                 229
# 6:        D01    P01-A01            LIG39   P01                 410
 
# > dim(merged)
# [1] 268073     15

merged <- merged %>% mutate_at(vars("hash_umis"), ~replace_na(., 0)) # add 0 if a cell has no hash
# > dim(merged)
# [1] 268073     15

#### Split top oligo by "." into top oligo 1, top oligo2...
#### Not actually accurate since it labelled as if these are ranked by top oligos 
#### This is actually meta data such as age, drug, etc 


suppressMessages(merged <- splitstackshape::cSplit(merged, "top_oligo", "."))
#                      cell hash_umis          pval          qval top_to_second_best_ratio Size_Factor n.umi perc_mitochondrial_umis P5_barcode P7_barcode RT_barcode
# 1: A05_D01_P01-A01_LIG107        11  5.218253e-29  1.002896e-24                 1.718059    37.46366   118               0.8474576        A05        D01    P01-A01
# 2: A05_D01_P01-A01_LIG160        38 6.462316e-289 1.195186e-283                11.537028   504.80688  1590               5.4716981        A05        D01    P01-A01
# 3: A05_D01_P01-A01_LIG163        20  0.000000e+00  0.000000e+00                 3.302363   125.72549   396               3.2828283        A05        D01    P01-A01
# 4:  A05_D01_P01-A01_LIG19         7  1.477237e-10  7.544251e-07                 2.496314    59.05288   186               0.0000000        A05        D01    P01-A01
# 5:  A05_D01_P01-A01_LIG23        11  2.701848e-82  2.056836e-77                 1.213414   103.81877   327               0.0000000        A05        D01    P01-A01
# 6:  A05_D01_P01-A01_LIG39        35 3.664147e-180 5.365227e-175                 1.741385   150.80709   475               0.4210526        A05        D01    P01-A01
# Ligation_barcode plate num_genes_expressed top_oligo_1 top_oligo_2 top_oligo_3 top_oligo_4 top_oligo_5 top_oligo_6
# 1:           LIG107   P01                 101       ATOH7          KO          72         hpf         P02          D8
# 2:           LIG160   P01                 878          WT          72         hpf         P02          C3        <NA>
#   3:           LIG163   P01                 274      WntC59          72          28           7          P1          A7
# 4:            LIG19   P01                 151       ATOH7          KO          72         hpf         P02          D4
# 5:            LIG23   P01                 229        DMSO          72          28           5          P1          G5
# 6:            LIG39   P01                 410       TGFB2          72          28           0          P1         F11

# > dim(merged)
# [1] 268073     20

### merge meta data from "merged" data frame to cds 
### can probably do this in a way where a separate "merged" data frame doesn't need to be created 
### Looks like they aren't using the p and q values or Sizefactor 
# this can be adapted based on how the hash_name is organized
### Need to figure out a standardized way to name hashes????
### Or not include this step and be done manually 
### or meta data can be provided 
colData(cds)$top_to_second_best_ratio = merged$top_to_second_best_ratio
colData(cds)$hash_umis = merged$hash_umis
colData(cds)$drug = merged$top_oligo_1
colData(cds)$age = merged$top_oligo_2
colData(cds)$temp = merged$top_oligo_3
colData(cds)$pheno = merged$top_oligo_4
colData(cds)$hash_plate = merged$top_oligo_5
colData(cds)$hash_well = merged$top_oligo_6

# drop any cells with no good hash assignment

#### Create an option where the threshold can be changed 
cds_filt <- cds[,colData(cds)$hash_umis > args$hash_umi_cutoff]
cds_filt <- cds_filt[,colData(cds_filt)$top_to_second_best_ratio > args$hash_ratio]

#
#                                          cell Size_Factor     n.umi perc_mitochondrial_umis   pcr_plate rt_plate_well    lig_well top_to_second_best_ratio
#                                   <character>   <numeric> <numeric>               <numeric> <character>   <character> <character>                <numeric>
# A05_D01_P01-A01_LIG160 A05_D01_P01-A01_LIG160    504.8069      1590                5.471698         D05       P01-A01      LIG160                 11.53703
# A05_D01_P01-A01_LIG163 A05_D01_P01-A01_LIG163    125.7255       396                3.282828         D05       P01-A01      LIG163                  3.30236
# A05_D01_P01-A01_LIG52   A05_D01_P01-A01_LIG52     34.6063       109                0.917431         D05       P01-A01       LIG52                      Inf
# A05_D01_P01-A01_LIG6     A05_D01_P01-A01_LIG6    208.5900       657               14.611872         D05       P01-A01        LIG6                  2.87628
# A05_D01_P01-A01_LIG66   A05_D01_P01-A01_LIG66    233.0366       734                5.585831         D05       P01-A01       LIG66                  3.05947
# A05_D01_P01-A02_LIG121 A05_D01_P01-A02_LIG121     97.4690       307                0.325733         D05       P01-A02      LIG121                  2.62458
#                        hash_umis        drug         age        temp       pheno  hash_plate   hash_well
#                        <numeric> <character> <character> <character> <character> <character> <character>
# A05_D01_P01-A01_LIG160        38          WT          72         hpf         P02          C3          NA
# A05_D01_P01-A01_LIG163        20      WntC59          72          28           7          P1          A7
# A05_D01_P01-A01_LIG52          7          WT          72         hpf         P02          C6          NA
# A05_D01_P01-A01_LIG6          31       ATOH7          KO          54         hpf         P02          B8
# A05_D01_P01-A01_LIG66         20       ATOH7          KO          54         hpf         P02          B4
# A05_D01_P01-A02_LIG121        24       ATOH7          KO          54         hpf         P02          B5

saveRDS(cds_filt,file=paste0(sample_name, "_hash_cds.RDS"))
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(garnett)
})

parser = argparse::ArgumentParser(description='Script to apply garnett.')
parser$add_argument('cds_path', help='File with cds.')
parser$add_argument('garnett_path', help='Csv with sample name and path to garnett model.')
parser$add_argument('sample_id', help='Sample ID.')
args = parser$parse_args()

garnett_file <- read.csv(args$garnett_path, header=FALSE)
garnett_file$orig <- garnett_file$V1
garnett_file$V1 <- gsub("-", "\\.", garnett_file$V1)
garnett_file$V1 <- gsub("_", "\\.", garnett_file$V1)
garnett_file$V1 <- gsub(" ", "\\.", garnett_file$V1)
garnett_file$V1 <- gsub("/", "\\.", garnett_file$V1)

if (sample_id %in% garnett_file$orig) {
    classifier_path <- garnett_file[garnett_file$orig == sample_id, ]$V2[1]
} else if (sample_id %in% garnett_file$V1) {
    classifier_path <- garnett_file[garnett_file$V1 == sample_id, ]$V2[1]
} else {
    classifier_path <- "NONE"
}

if (classifier_path == "NONE") {
    file.copy(args$cds_path, "new_cds/")
} else {
    cds <- readRDS(args$cds_path)
    classifier <- readRDS(classifier_path)

    cds <- classify_cells(cds, classifier,
                           db = "none",
                           cluster_extend = FALSE,
                           cds_gene_id_type = "ENSEMBL")
    
    saveRDS(cds, file=paste0("new_cds/", sample_id, "_cds.RDS"))

}
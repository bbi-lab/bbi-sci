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

if (args$sample_id %in% garnett_file$orig) {
    classifier_path <- garnett_file[garnett_file$orig == args$sample_id, ]$V2
} else if (args$sample_id %in% garnett_file$V1) {
    classifier_path <- garnett_file[garnett_file$V1 == args$sample_id, ]$V2
} else {
    classifier_path <- "NONE"
}

fileConn<-file("garnett_error.txt")
if (classifier_path == "NONE") {
    file.copy(args$cds_path, "new_cds/")
} else {
    tryCatch({
    cds <- readRDS(args$cds_path)
    for (val in classifier_path) {
        classifier <- readRDS(as.character(val))
        classifier_name <- unlist(stringr::str_split(val, "/"))
        classifier_name <- classifier_name[length(classifier_name)]
        classifier_name <- gsub(".RDS", "", classifier_name)

        cds <- classify_cells(cds, classifier,
                            db = "none",
                            cluster_extend = FALSE,
                            cds_gene_id_type = "ENSEMBL")
        
        names(colData(cds))[names(colData(cds)) == "cell_type"] <- paste0("garnett_type_", classifier_name)
    }
    saveRDS(cds, file=paste0("new_cds/", args$sample_id, "_cds.RDS"))
    writeLines("ok", fileConn)
    }, error = function(e) {
        file.copy(args$cds_path, "new_cds/")
        writeLines(as.character(e), fileConn)
    })
 close(fileConn)
}

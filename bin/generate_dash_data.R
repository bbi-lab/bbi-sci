#!/usr/bin/env Rscript

library(jsonlite)
library(argparse)

parser = argparse::ArgumentParser(description='Script to generate dashboard data file.')
parser$add_argument('dup_file', help='File of all duplication stats.')
parser$add_argument('output_folder', help='Output folder.')
parser$add_argument('cell_counts', help='Cell counts file.')
parser$add_argument('barn_col', help='File of concatenated collision values.')
parser$add_argument('garnett_csv', help='File of Garnett models or false.')
args = parser$parse_args()

all_dups <- read.csv(args$dup_file, header=FALSE, stringsAsFactors=FALSE)
output_folder <- args$output_folder
count_info <- read.table(args$cell_counts, stringsAsFactors=FALSE)

project_name <- unlist(stringr::str_split(output_folder, "/"))
project_name <- project_name[[length(project_name)]]

c100 <- sum(count_info[count_info$V2 == 100,]$V3)
c500 <- sum(count_info[count_info$V2 == 500,]$V3)
c1000 <- sum(count_info[count_info$V2 == 1000,]$V3)

count_info_tab <- count_info

count_info_tab$V1 <- gsub("_cell_qc.csv", "", count_info_tab$V1)
ct100 <- count_info_tab[count_info_tab$V2 == 100,]
ct500 <- count_info_tab[count_info_tab$V2 == 500,]
ct1000 <- count_info_tab[count_info_tab$V2 == 1000,]
row.names(ct100) <- ct100$V1
row.names(ct500) <- ct500$V1
row.names(ct1000) <- ct1000$V1

all_dups$V1 <- as.character(all_dups$V1)
all_dups$c100 <- ct100[all_dups$V1,"V3"]
all_dups$c1000 <- ct1000[all_dups$V1,"V3"]

all_dups$V5[all_dups$V7 > 0] <- "Fail"
all_dups$V6[all_dups$V7 > 0] <-  "Fail"
all_dups$V6[all_dups$V6 == "NaN%"] <-  "Fail"

row.names(all_dups) <- all_dups$V1
names(all_dups) <- c("Sample", "Total_reads",
                     "Total_UMIs",
                     "Duplication_rate",
                     "Doublet_Number", 
                     "Doublet_Percent",
                     "Doublet_NAs",
                     "Cells_100_UMIs",
                     "Cells_1000_UMIs" 
                     )

all_dups$Doublet_Number[is.na(all_dups$Doublet_Number)] <- "Fail"
all_dups$Doublet_Percent[is.na(all_dups$Doublet_Percent)] <- "Fail"

if (args$garnett_csv != "false") {
  garnett_file <- read.csv(args$garnett_csv, header=FALSE, stringsAsFactors = FALSE)
  garnett_file$V2 <- gsub(".RDS", "", garnett_file$V2)
  garnett_file$V2 <- sapply(garnett_file$V2, function(x) {
     y <- unlist(stringr::str_split(x, "/"))
     y[length(y)]
  })
  all_dups$Garnett_model <- NA
  garnett_file$samp_fixed <- garnett_file$V1
  garnett_file$samp_fixed <- gsub("_", ".", garnett_file$samp_fixed)
  garnett_file$samp_fixed <- gsub("-", ".", garnett_file$samp_fixed)
  garnett_file$samp_fixed <- gsub(" ", ".", garnett_file$samp_fixed)
  for (samp in all_dups$Sample) {
    if (samp %in% garnett_file$V1 | samp %in% garnett_file$samp_fixed) {
      all_dups$Garnett_model[all_dups$Sample == samp] <- list(garnett_file[garnett_file$V1 == samp | garnett_file$samp_fixed == samp, "V2"])
    
      if (all_dups$Cells_100_UMIs[all_dups$Sample == samp] == 0) {
        all_dups$Garnett_model[all_dups$Sample == samp] <- "no_cells"
      }
    }
  }
}


all_dup_lst <- apply(all_dups, 1, as.list)
sample_list <- as.character(all_dups$Sample)[order(as.character(all_dups$Sample))]

if("Sentinel" %in% sample_list) {
  sample_list <- c("Sentinel", setdiff(sample_list, c("Sentinel")))
}
barn_collision <- NA
if("Barnyard" %in% sample_list) {
  sample_list <- c("Barnyard", setdiff(sample_list, c("Barnyard")))
  barn_collision <- read.table(args$barn_col)
  barn_collision <- barn_collision[barn_collision$V1 == "Barnyard", "V2"]
}
sample_list <- as.list(sample_list)
json_info <- list("run_name" = project_name,
                  "cell_counts" = c(sum(all_dups$Cells_100_UMIs), sum(all_dups$Cells_1000_UMIs)),
                  "sample_list" = sample_list,
                  "barn_collision" = barn_collision,
                  "sample_stats" = all_dup_lst)
                  
fileConn<-file("data.js")
writeLines(c("const run_data =", toJSON(json_info, na='null',  pretty=TRUE, auto_unbox=TRUE)), fileConn)
close(fileConn)

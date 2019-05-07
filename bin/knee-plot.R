#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(methods)
    library(dplyr)
    library(ggplot2)
    library(mclust)
    library(argparse)
})

parser = argparse::ArgumentParser(description='Script to do cell calling and output a threshold for cell calls.')
parser$add_argument('input_file', help='File with per cell exon + intron counts.')
parser$add_argument('--specify_cutoff', type='integer', help='Optional. Specifies a cutoff rather than choosing a UMI cutoff automatically.')
parser$add_argument('--umi_count_threshold_file', required=TRUE, help='Output file that cutoff value is written to.')
parser$add_argument('--knee_plot', required=TRUE, help='File to save knee plot to (used to show UMI per cell cutoff).')
args = parser$parse_args()

cutoff = NULL
if (! is.null(args$cutoff)) {
    cutoff = args$cutoff
}

plot_name <- gsub("%", "%%", args$knee_plot)

df = read.table(
    args$input_file,
    col.names = c("sample", "barcode", "n.umi"),
    colClasses = c("factor", "character", "integer"))

# if empty return empty plot
if (nrow(df) == 0) {
df <- data.frame()
empty <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
ggsave(plot_name, plot = empty, units = "in", width = 3.5*1.3, height = 3.5)

fileConn<-file(args$umi_count_threshold_file)
writeLines(c(as.character(0)), fileConn)
close(fileConn)
} else {
# Get the two populations with mclust
if(is.null(cutoff)) {
    background_call = Mclust(data.frame(log10(df$n.umi)),G=c(2))
    cutoff = min(df[which(background_call$classification == 2 & background_call$uncertainty < 0.005),3])
}

# Output plot
df = df %>% mutate(n.umi.rank = min_rank(-n.umi))

plot = ggplot(df %>%
                arrange(-n.umi) %>%
                select(n.umi, n.umi.rank)
                %>% distinct(),
            aes(x = n.umi.rank, y = n.umi)) +
        geom_line(size = 0.8) +
        scale_x_log10() +
        scale_y_log10() +
        xlab("# of barcodes") +
        ylab("UMI count threshold") +
        theme_bw()

    plot = plot +
        geom_hline(yintercept = cutoff, size = 1.2, color = "firebrick2")

ggsave(plot_name, plot = plot, units = "in", width = 3.5*1.3, height = 3.5)

# Write out threshold to file
fileConn<-file(args$umi_count_threshold_file)
writeLines(c(as.character(cutoff)), fileConn)
close(fileConn)
}

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(methods)
    library(dplyr)
    library(ggplot2)
})


parser = argparse::ArgumentParser(description='Script to generate knee plot for hash umis.')
parser$add_argument('hash_umi', help='File with cds.')
parser$add_argument('sample_name', help='Sample name')
args = parser$parse_args()

df = read.table(
    args$hash_umi,
    col.names = c("sample", "barcode", "n.umi"),
    colClasses = c("factor", "character", "integer"))

# df = data.frame(args$hash_umi)
# colnames(df) <- c("sample", "barcode", "n.umi")
# colClasses(df) <- c("factor", "character", "integer")

df = df %>% group_by(sample) %>% mutate(n.umi.rank = min_rank(-n.umi)) %>% ungroup()


if (dim(df)[1] != 0) {   
    for (this.sample in levels(df$sample)) {
        plot = ggplot(
            df %>% filter(sample == this.sample) %>%
                arrange(-n.umi) %>% select(n.umi, n.umi.rank) %>% distinct(),
            aes(x = n.umi.rank, y = n.umi)) +
            geom_line(size = 0.8) +
            # scale_x_log10(limits = c(10, NA),
            #               breaks = c(10, 50, 100, 500, 1000,2000,4000, 8000, 16000, 32000, 64000, 128000,250000,500000,750000, 1000000)) +
            # scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000)) +
            scale_x_log10() +
            scale_y_log10() +
            xlab("# of barcodes") +
            ylab("hash UMI count threshold") +
            theme_bw()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

        # if (length(args) >= 3) {
        #     plot = plot +
        #         geom_hline(yintercept = cutoff, size = 1.2, color = "firebrick2")
        # }

            # plot = plot +x
            #     geom_hline(yintercept = cutoff, size = 1.2, color = "firebrick2")
        
        # ggsave(paste(args[2], "/", this.sample, ".pdf", sep = ""),
        #     plot = plot, units = "in", width = 3.5*1.618, height = 3.5)


        ggsave(paste0(args$sample_name, "_hash_knee_plot.png"),
            plot = plot, units = "in", width = 3.5*1.618, height = 3.5)

        # cat(this.sample, "\n", sep ="")
    }
} else {

    file_name <- paste0(args$sample_name, "_hash_knee_plot.png")
    ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient Hash UMIs")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="RT Barcodes", y = "hash UMI count threshold")
    ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')
}



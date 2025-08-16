#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(viridis)
    library(data.table)
    library(dplyr)
    library(BPCells)
})


parser = argparse::ArgumentParser(description='Script to pseudobulk plots')
parser$add_argument('cds_path', help='File with cds.')
parser$add_argument('sample_name', help='Sample name.')
args = parser$parse_args()

# Get average gene expression for each sample and RT barcode 
calc_pseudobulk <- function(cds) {
  psb_list <- list()
    
  # Iterate through each RT well and get average gene expression across RTs for each gene
  for (rt in levels(as.factor(as.character(cds$RT_barcode)))) {
    avg <- unname(rowMeans(normalized_counts(cds)[,cds$RT_barcode==rt]))
    psb_list[[rt]] <- avg
  }
  
  pseudo_bulkDF <- data.frame(psb_list)
  rownames(pseudo_bulkDF) <- paste(rowData(cds)$gene_short_name, rowData(cds)$id, sep="_") ## Set gene names as rownames 
  
  ## R correlation using Pearson correlation and R^2
  cor_DF <- cor(pseudo_bulkDF)^2
  psb <- melt(cor_DF, id.vars=c(rownames(cor_DF)))
  colnames(psb) <- c("rt1", "rt2", "pearson_correlation")
  
  write.table(psb, paste0(args$sample_name, "_pseudobulk_rt_correlations.txt"), sep="\t", quote=FALSE, row.names = FALSE)
  
  return(psb)
}

tryCatch({
    cds <- load_monocle_objects(args$cds_path)

    psb <- calc_pseudobulk(cds)
    summary(psb$pearson_correlation)

    ### Create heatmap for all samples
    heatmap <- ggplot(psb, aes(x=rt1, y=rt2, fill=pearson_correlation)) +
    geom_tile() +
    scale_fill_viridis(direction=-1, limits=c(0.60,1.000)) +
    theme_minimal(base_size=14) + 
    theme(axis.text.x=element_text(angle=90, hjust=1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    # scale_fill_continuous(low = "yellow", high = "purple", limits = c(0.5, 1.000))

    ggsave(paste0(args$sample_name, "_correlation_heatmap.png"), heatmap)

    hist_plot <- ggplot(psb, aes(x = pearson_correlation)) +
    geom_histogram(binwidth = .005, fill = "grey", color = "black") +
    theme_bw(base_size=20) +
    # coord_cartesian(xlim = c(min(psb$value), 1)) +
    coord_cartesian(xlim = c(.60, 1)) +
    geom_vline(xintercept = 0.90, color = "red", linetype = "dashed", linewidth = 0.5) 

    ggsave(paste0(args$sample_name, "_correlation_hist.png"), heatmap)

}, error = function(e) {
    
    ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for RT correlation")) + theme(legend.position = "none")
    ggsave(paste0(args$sample_name, "_correlation_hist.png"), ggp_obj, device='png', width=5, height=5, dpi=600, units='in')
    ggsave(paste0(args$sample_name, "_correlation_heatmap.png"), ggp_obj, device='png', width=5, height=5, dpi=600, units='in')
    write.table("Insufficent Cells for RT correlation", file=paste0(args$sample_name, "_pseudobulk_rt_correlations.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

})
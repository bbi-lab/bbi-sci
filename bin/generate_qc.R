#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(methods)
  library(dplyr)
  library(ggplot2)
  library(mclust)
  library(argparse)
  library(monocle3)
  library(tidyverse)
  library(cowplot)
  library(randomcoloR)
  library(scales)
  library(reshape2)
  library(viridis)
})

parser = argparse::ArgumentParser(description='Script to generate qc plots.')
parser$add_argument('cds_path', help='File with cds.')
parser$add_argument('umis_file', help='File with umis_per_barcode.')
parser$add_argument('sample_name', help='Sample name.')
parser$add_argument('empty_drops', help='RDS file from emptyDrops.')
parser$add_argument('hash', help='hash run or not.')
parser$add_argument('--specify_cutoff', type='integer', help='Optional. Specifies a cutoff rather than choosing a UMI cutoff automatically.')
args = parser$parse_args()

sample_name <- args$sample_name

# Read empty drops RDS object
emptydrops_data <- readRDS(args$empty_drops)


#################################################
##   FUNCTIONS FOR EXTRACTING BARCODES AND     ##
##          PERFORMING QC CHECKS               ##
#################################################



# Extract column data from monocle object and separate 
# barcodes and plate info from the cell name 
# Returns a new cds object with barcode and plate info
extractBarcode <- function(cds) {
    df <- as.data.frame(colData(cds))

    meta_types <- c("P5_barcode", "P7_barcode", "RT_barcode", "Ligation_barcode")
    meta <- separate(df, cell, into=meta_types, sep="_", remove=FALSE)

    # Fill in barcode information in cds object from meta object 
    for (m in meta_types) {
        colData(cds)[,m] <- meta[[m]]
    }

    # Extract plate number
    cds$RT_plate <- sapply(strsplit(as.character(meta$RT_barcode), "-"), `[`, 1)

    lig_info <- list(P01=1:96,
                 P02=97:192,
                 P03=193:288,
                 P04=289:384)

    # Convert barcode to numeric
    lig_nums <- as.integer(gsub("LIG", "", cds$Ligation_barcode))

    # Map each ligation barcode to its plate
    cds$Ligation_plate <- vapply(lig_nums, function(n) {
      name <- names(lig_info)[sapply(lig_info, function(x) n %in% x)]
      if (length(name) == 0) return(NA_character_) else return(name)
    }, character(1))

    return (cds)
}

# Function to return labels/text for stat_summary
# Returns labels for ymax data (here, it is max number of cells)
n_fun <- function(y){return(data.frame(y=max(y), label = paste0(length(y))))}

##############################
####   UMI Violin plots   ####
##############################

# Generate and returns a combined violin plot for umis by plates and p5 

generate_umi_plots <- function(cds) {
  n_fun <- function(y){return(data.frame(y=max(y), label = paste0(length(y))))}

  barcode_list <- list("RT_plate", "Ligation_plate", "P5_barcode")
  umi_plot_list <- list()

  plot_list <- lapply(barcode_list, function(barcode) {

    barcode_name <- gsub("_", " ", barcode)

    # Plot UMIs by RT, Ligation and P5 plates 
    umi_plot <- ggplot(data.frame(colData(cds)), aes(x=.data[[barcode]], y=n.umi)) +
      geom_violin(fill="salmon")+
      geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
      theme_light(base_size=10) +
      scale_y_log10(labels=trans_format("log10", math_format(10^.x)),
                    limits=c(100, max(cds$n.umi) *30 )) +
      xlab(barcode_name) +
      ylab("UMIs") +
      stat_summary(fun.data = n_fun, geom = "text", angle=90, hjust = -0.1, size=3) +
      theme(legend.position="none") +
      theme(axis.text.x = element_text(angle=90, hjust=1),
            axis.title.x = element_text(margin=margin(t=10)),
            axis.title.y= element_text(margin=margin(r=10))) 

      umi_plot_list[[barcode]] <- umi_plot

    return(umi_plot_list)
    })

  flat_plot_list <- unlist(plot_list, recursive=FALSE)

  # Add spacing between combined ggplots
  plot_spacer <- function() ggplot() + theme_void()
  combined_umi <- plot_grid(flat_plot_list$RT_plate,
                            plot_spacer(),
                            flat_plot_list$Ligation_plate,
                            plot_spacer(), 
                            flat_plot_list$P5_barcode,
                            align="hv", 
                            ncol=1,
                            rel_heights = c(1, 0.25, 1, 0.25, 1))

  return(combined_umi) 
}


##############################
####   Pseudobulk plots   ####
##############################

# Calculate pseudobulk gene expression correlations for each sample 
# by barcodes and plates
# Outputs a combined correlation heatmap, combined distribution histogram, 
# and writes out tables for each barcode/plate correlations

calc_pseudobulk <- function(cds, sample_name) {
  
  barcode_list <- list("RT_barcode", "Ligation_plate", "P5_barcode")
  
  norm_counts <- normalized_counts(cds)   
  plot_list = list()
  
  plot_list <- lapply(barcode_list, function(barcode) {
    print(barcode)
    
    psb_list <- list()
    
    # Precompute barcode groupings
    barcode_vec <- as.character(cds[[barcode]])
    well_groups <- split(seq_along(barcode_vec), barcode_vec)
    
    # Restrict to first 10 wells only
    wells_to_use <- names(well_groups)
    
    # Use lapply with Matrix-aware rowMeans
    psb_list <- lapply(wells_to_use, function(well) {
      cols <- well_groups[[well]]
      Matrix::rowMeans(norm_counts[, cols, drop = FALSE])
    })
    
    names(psb_list) <- wells_to_use
    
    pseudo_bulkDF <- data.frame(psb_list)
    rownames(pseudo_bulkDF) <- paste(rowData(cds)$gene_short_name, rowData(cds)$id, sep="_") ## Set gene names as rownames
    
    ## R correlation using Pearson correlation and R^2
    cor_DF <- cor(pseudo_bulkDF)^2

    psb <- melt(cor_DF, id.vars=c(rownames(cor_DF)))

    barcode_name <- gsub("_", " ", barcode)
    colnames(psb) <- c(paste0(barcode_name," 1"), paste0(barcode_name, " 2"), "pearson correlation")

    ### Create heatmap for all samples
    
    heatmap_theme <- NULL 
    
    if (barcode=="Ligation_plate") {
      heatmap_theme <- theme(axis.text.x = element_blank(), 
                             axis.text.y = element_blank(),
                             legend.text = element_text(size = 6),   
                             legend.title = element_text(size = 6)) 
    } else {
      heatmap_theme <- theme(axis.ticks = element_line(linewidth = 0.2),
                              axis.ticks.length = unit(0.1, "cm"),
                              axis.title.x = element_text(size = 8, 
                                                          margin=margin(t=6)),
                              axis.title.y = element_text(size = 8,
                                                          margin=margin(r=10)),
                              axis.text.x = element_text(size=6, angle=90),
                              axis.text.y = element_text(size=6), 
                              legend.text = element_text(size = 6),   
                              legend.title = element_text(size = 6)) 
    }
    heatmap <- ggplot(psb, aes(x=psb[,1], y=psb[,2], fill=.data[["pearson correlation"]])) +
      geom_tile() +
      coord_fixed() +
      scale_fill_viridis(direction=-1, limits=c(0,1.000)) +
      theme_minimal(base_size=8) +
      xlab(colnames(psb)[1]) + 
      ylab(colnames(psb)[2]) +
      heatmap_theme 
  
    plot_list[[paste0("heatmap_",barcode)]] <- heatmap
    
    
    ### Create distribution of barcode correlations
    hist_plot <- ggplot(psb, aes(x = .data[["pearson correlation"]])) +
      geom_histogram(binwidth = .003, fill = "grey", color = "black") +
      theme_bw(base_size=12) +
      # coord_cartesian(xlim = c(min(psb$value), 1)) +
      coord_cartesian(xlim = c(.60, 1)) +
      geom_vline(xintercept = 0.90, color = "red", linetype = "dashed", linewidth = 0.5) + 
      ggtitle(barcode_name)
    
    
    plot_list[[paste0("hist_",barcode)]] <- hist_plot

    write.table(psb, paste0(sample_name, "_", barcode, "_pseudobulk_correlations.txt"), sep="\t", quote=FALSE, row.names = FALSE)
    
    return(plot_list)
    
  })
  
  
  # Flatten one level so you get a named list of ggplot objects

  flat_plot_list <- unlist(plot_list, recursive = FALSE)
  flat_plot_list
  
  plot_spacer <- function() ggplot() + theme_void()
  
  combined_heatmap <- plot_grid(flat_plot_list$heatmap_RT_barcode, 
                              plot_spacer(),
                              flat_plot_list$heatmap_Ligation_plate,
                              plot_spacer(),
                              flat_plot_list$heatmap_P5_barcode,
                              ncol=1, 
                              rel_heights = c(1, 0.25, 1, 0.25, 1),
                              align="hv")
  
  
  combined_hist <- plot_grid(flat_plot_list$hist_RT_barcode, 
                                 plot_spacer(),
                                 flat_plot_list$hist_Ligation_plate,
                                 plot_spacer(), 
                                 flat_plot_list$hist_P5_barcode,
                                 ncol=1, 
                                 rel_heights = c(1, 0.25, 1, 0.25, 1),
                                 align="hv")
  
  combined_hist
  
  plots <- list(combined_heatmap, combined_hist)

  return (plots)
  
}


##############################
####    Genes by UMIs     ####
##############################

# Generate and return a combined genes by umi plot colored by emptyDrops
# and by percent mitochondrial umis

generate_genes_by_umis <- function(cds) {

  ### By perc mito

  genes_by_umi1 <- ggplot(data.frame(colData(cds)), aes(x=n.umi, y=num_genes_expressed, color=perc_mitochondrial_umis)) +
    geom_point(aes(alpha=0.3), size=0.5) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          aspect.ratio = 1,
          legend.title = element_text(size = 8),
          legend.position="right") +
    xlab("UMIs") +
    ylab("Number of Genes Captured") +
    scale_y_log10(labels=trans_format("log10", math_format(10^.x)), limit=c(1,max(cds$num_genes_expressed) + 1000)) +
    scale_x_log10(labels=trans_format("log10", math_format(10^.x))) +
    geom_abline(slope=1, color="grey") +
    scale_colour_continuous(limit=c(0, 100), low = "black", high = "green") + 
    guides(alpha = "none", size="none")

  ### By emptydrops 

  genes_by_umi2 <- ggplot(data.frame(colData(cds)), aes(x=n.umi, y=num_genes_expressed, color=emptyDrops_FDR)) +
    geom_point(aes(alpha=0.3), size =0.5) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          aspect.ratio = 1,
          legend.title = element_text(size = 8),
          legend.position="right") +
    xlab("UMIs") +
    ylab("Number of Genes Captured") +
    scale_y_log10(labels=trans_format("log10", math_format(10^.x)), limit=c(1,max(cds$num_genes_expressed) + 1000)) +
    scale_x_log10(labels=trans_format("log10", math_format(10^.x))) +
    geom_abline(slope=1, color="grey") + 
    guides(alpha = "none", size = "none")

  genes_by_umi <- plot_grid(genes_by_umi1, genes_by_umi2, align="hv", ncol=1)
  ggsave("genes_by_umi.png", genes_by_umi, width=5, height=7, dpi=600, units='in')

  return(genes_by_umi)

}


##############################
####    Hash Knee Plots   ####
##############################


# Generate ranked barcodes by hash umis knee plot 
# and hash top-second-best-ratio knee plot
# Returns a combined plot

generate_hash_plots <- function(cds) {

  print("here in hash")
  df <- data.frame(colData(cds))

  # Rank barcodes by hash umis
  df = df %>% mutate(n.umi.rank = min_rank(-hash_umis)) %>% ungroup() %>% data.frame()

  # Hash umi by barcode rank knee plot 
  # Removes duplicate ranks
  hash_umis <- ggplot(
    df %>% 
      arrange(-hash_umis) %>% select(hash_umis, n.umi.rank) %>% distinct(),
    aes(x = n.umi.rank, y = hash_umis)) +
    geom_point(shape = 21) +
    scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + 
    scale_y_log10(labels=trans_format("log10", math_format(10^.x))) + 
    xlab("# of barcodes") +
    ylab("hash umis") +
    theme_bw(base_size = 10)+
    theme(axis.title.y = element_text(margin = margin(r = 15))) 


  # Hash top-to-second-best ratio vs hash barcode ranks
  df <- data.frame(colData(cds))

  df <- df[!is.na(df$hash_umis) &
            df$hash_umis > 0 &
            is.finite(df$top_to_second_best_ratio),]

  df = df %>% mutate(n.umi.rank = min_rank(-top_to_second_best_ratio)) %>% ungroup() %>% data.frame()
  head(df)

  # Hash top-to-second best ratio vs hash barcode ranks knee plot

  hash_top_knee <- ggplot(
    df %>% 
      arrange(-top_to_second_best_ratio) %>% select(top_to_second_best_ratio, n.umi.rank) %>% distinct(),
    aes(x = n.umi.rank, y = top_to_second_best_ratio)) +
    geom_point(shape=21) +
    scale_x_log10(labels=trans_format("log10", math_format(10^.x))) + 
    scale_y_log10(labels=trans_format("log10", math_format(10^.x))) + 
    xlab("# of barcodes") +
    ylab("top to second best ratio") +
    theme_bw(base_size = 10)+
    theme(axis.title.y = element_text(margin = margin(r = 15))) + 
    geom_hline(yintercept = 2, linetype= "dotted", color = "firebrick", size= 0.5) + 
    geom_hline(yintercept = 5, linetype= "dotted", color = "cornflowerblue", size =0.5) +
    annotate("text", x = 100 , y = 2 * 0.45, label = "top-to-second-best cutoff = 2", vjust = -0.5, hjust = 1, size = 2) +
    annotate("text", x = 100  , y = 5 * 1.10, label = "top-to-second-best cutoff = 5", vjust = -0.5, hjust = 1, size = 2) 
    

  hash_top_knee


  combined_knee <- plot_grid(hash_umis, hash_top_knee, ncol=1, align = "hv")

  return(combined_knee)
}


#################################################
##             START QC CHECKS                 ##
#################################################


gen_plots <- function(sample_name, sample_path) {
  
  # if (file.size(paste0(sample_path,"/bpcells_matrix_dir/col_names")) != 0) {
  #   samp_cds <- load_monocle_objects(sample_path)
  # } else {
  #   samp_cds <- readRDS(paste0(sample_path,"/cds_object.rds"))
  # }

  # garnett_mods <- names(colData(samp_cds))[grepl("garnett_type", names(colData(samp_cds)))]  

  samp_cds <- tryCatch({

    samp_cds <- load_monocle_objects(sample_path)

    # If Empty drops was ran, then add info for Empty Drops FDR <= 0.01
    if(is(emptydrops_data, 'DFrame')) {
      samp_cds$emptyDrops_FDR_0.01 <- ifelse(
        is.na(colData(samp_cds)$emptyDrops_FDR), "NA", 
        ifelse(colData(samp_cds)$emptyDrops_FDR <= 0.01, "TRUE", "FALSE"))
    }

    samp_cds <- extractBarcode(samp_cds)
    samp_cds <- detect_genes(samp_cds)
    
    umi_plots <- generate_umi_plots(samp_cds)
    ggsave(paste0(sample_name, "_umi.png"), umi_plots, width=5, height=8, dpi=600, units="in")

    genes_umi_plot <- generate_genes_by_umis(samp_cds)
    ggsave(paste0(sample_name, "_genes_by_umi.png"), genes_umi_plot, width=5, height=7, dpi=600, units='in')

    if(is(emptydrops_data, 'DFrame')) {
      keep_na <- samp_cds[,is.na(colData(samp_cds)$emptyDrops_FDR)]
      keep_cells <- samp_cds[,!is.na(colData(samp_cds)$emptyDrops_FDR) & colData(samp_cds)$emptyDrops_FDR <= 0.01]
      samp_cds <- combine_cds(list(keep_na, keep_cells), sample_col_name="og_cds")
    }

    psb_plots <- calc_pseudobulk(samp_cds, sample_name)
    ggsave(paste0(sample_name, "_pseudobulk_heatmap.png"), psb_plots[[1]], width=5, height=10, dpi=800, units='in')
    ggsave(paste0(sample_name, "_pseudobulk_histogram.png"), psb_plots[[2]], width=5, height=10, dpi=800, units='in')

    # Skip hash plots if not a hash run 
    if (args$hash != 'false') {
      hash_plots <- generate_hash_plots(samp_cds)
      ggsave(paste0(sample_name, "_hash_plots.png"), hash_plots, width=5, height=7, dpi=600, units='in')

    } else {
      file_name <- paste0(sample_name, "_hash_plots.png")
      ggp_obj <- ggplot() + theme_void() + geom_text(aes(x = 1, y = 1, label = "Process hash skipped.")) 
      ggsave(filename=file_name, ggp_obj, device='png')

    }

  }, error = function(e) {

    # Generate empty plots if error
    file_name <- paste0(sample_name, "_umi.png")
    ggp_obj <- ggplot() + theme_void() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for UMI plot")) 
    ggsave(filename=file_name, ggp_obj, device='png')

    file_name <- paste0(sample_name, "_pseudobulk_heatmap.png")
    ggp_obj <- ggplot() + theme_void() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for Pseudobulk correlations")) 
    ggsave(filename=file_name, ggp_obj, device='png')

    file_name <- paste0(sample_name, "_pseudobulk_histogram.png")
    ggp_obj <- ggplot() + theme_void() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for Pseudobulk correlations")) 
    ggsave(filename=file_name, ggp_obj, device='png')

    if (args$hash != 'false') { 
      file_name <- paste0(sample_name, "_hash_plots.png")
      ggp_obj <- ggplot() + theme_void() + geom_text(aes(x = 1, y = 1, label = "Insufficient hash umis")) 
      ggsave(filename=file_name, ggp_obj, device='png')

    } else {
      file_name <- paste0(sample_name, "_hash_plots.png")
      ggp_obj <- ggplot() + theme_void() + geom_text(aes(x = 1, y = 1, label = "Process hash skipped.")) 
      ggsave(filename=file_name, ggp_obj, device='png')
    }


  })

}


# Generate UMAP, cell_qc
cds <- gen_plots(args$sample_name, args$cds_path)

# Knee plot

# cutoff = NULL
cutoff = args$specify_cutoff
gen_knee <- function(sample_name, cutoff) {

  
  df = read.table(
    args$umis_file,
    col.names = c("barcode", "n.umi"),
    colClasses = c("character", "integer"))
  
  # if empty return empty plot
  if (nrow(df) == 0) {
    df <- data.frame()
    empty <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
    ggsave(paste0(sample_name, "_knee_plot.png"), plot = empty, units = "in", width = 3.5*1.3, height = 3.5)
    
    fileConn<-file(paste0(sample_name, "_umi_cutoff.png"))
    writeLines(c(as.character(0)), fileConn)
    close(fileConn)
  } else {
    # Get the two populations with mclust
    # if(is.null(cutoff)) {
    #   background_call = Mclust(data.frame(log10(df$n.umi)),G=c(2))
    #   cutoff = min(df[which(background_call$classification == 2 & background_call$uncertainty < 0.005),3])
    # } else  {
    #   cutoff = args$cutoff
    # }

    if(is.null(cutoff)) {
      background_call = Mclust(data.frame(log10(df$n.umi)),G=c(2))
      cutoff = min(df[which(background_call$classification == 2 & background_call$uncertainty < 0.005),3])
    } 

    plot = NULL

    # Color barcode rank knee plot by empty drops FDR if empty drops was called 
    if(is(emptydrops_data, 'DFrame')) {

      df = cbind(df, emptyDrops_FDR=emptydrops_data$FDR)
      df$emptyDrops_FDR_0.01 <- ifelse(is.na(emptydrops_data$FDR), "NA", ifelse(df$emptyDrops_FDR <= 0.01, "TRUE", "FALSE"))
      df = df %>% mutate(n.umi.rank = min_rank(-n.umi))

      plot = ggplot(df %>%
                      arrange(-n.umi) %>%
                      select(n.umi, n.umi.rank, emptyDrops_FDR_0.01)
                    %>% distinct(),
                    aes(x = n.umi.rank, y = n.umi, color=emptyDrops_FDR_0.01)) +
        geom_point(size=1, shape=1, stroke=0.3) +
        scale_color_manual(values = c("TRUE" = "cornflowerblue", "FALSE" = "red", "NA"="grey")) +
        scale_x_log10() +
        scale_y_log10() +
        xlab("# of barcodes") +
        ylab("UMI count threshold") +
        theme_bw() + 
        theme(legend.position="bottom")
    } else {
    # Output plot
      df = df %>% mutate(n.umi.rank = min_rank(-n.umi))

      plot = ggplot(df %>%
                      arrange(-n.umi) %>%
                      select(n.umi, n.umi.rank)
                    %>% distinct(),
                    aes(x = n.umi.rank, y = n.umi)) +
        geom_point(size=1, shape=1, stroke=0.3, color = "black") +
        # geom_line(size = 0.8) +
        scale_x_log10() +
        scale_y_log10() +
        xlab("# of barcodes") +
        ylab("UMI count threshold") +
        theme_bw()
    }
    
    plot = plot +
      geom_hline(aes(yintercept = cutoff), linetype="dotted", linewidth = .5, color = "firebrick2") + 
      annotate("text", x = max(df$n.umi.rank), y = cutoff, label = "umi cutoff", vjust = -0.5, hjust = 1, size = 2) 
    
    ggsave(paste0(sample_name, "_knee_plot.png"), plot = plot, units = "in", width = 3.5*1.3, height = 3.5)
    
    # Write out threshold to file
    fileConn<-file(paste0(sample_name, "_umi_cutoff.txt"))
    writeLines(c(as.character(cutoff)), fileConn)
    close(fileConn)
  }
}


suppressMessages(gen_knee(args$sample_name, args$specify_cutoff))



# Generate Barnyard collision rates between mouse and human genes 

if (sample_name == "Barnyard") {

  fData(cds)$mouse <- grepl("ENSMUSG", fData(cds)$id)
  fData(cds)$human <- grepl("ENSG", fData(cds)$id)

  pData(cds)$mouse_reads <- Matrix::colSums(exprs(cds)[fData(cds)$mouse,])
  pData(cds)$human_reads <- Matrix::colSums(exprs(cds)[fData(cds)$human,])
  pData(cds)$total_reads <- pData(cds)$mouse_reads + pData(cds)$human_reads
  pData(cds)$human_perc <- pData(cds)$human_reads/pData(cds)$total_reads
  pData(cds)$mouse_perc <- pData(cds)$mouse_reads/pData(cds)$total_reads
  pData(cds)$collision <- ifelse(pData(cds)$human_perc >= .9 | pData(cds)$mouse_perc >= .9, FALSE, TRUE)


  plot = ggplot(as.data.frame(pData(cds)), aes(mouse_reads, human_reads, color = collision)) +
    geom_point(size = .8) +
    theme_bw() +
    scale_color_manual(values = c("black", "red")) +
    theme(legend.position = "none") +
    xlab("Mouse UMIs") +
    ylab("Human UMIs")

  ggsave("Barnyard_plot.png", plot = plot, units = "in", width = 3.5*1.3, height = 3.5)

  collision_rate <- round(sum(pData(cds)$collision/nrow(pData(cds))) * 200, 1)
  fileConn<-file("Barn_collision.txt")
  writeLines(paste0(args$sample_name, "\t", collision_rate, "%"), fileConn)
  close(fileConn)

} else {
    fileConn<-file(paste0(args$sample_name, "_no_collision.txt"))
    writeLines(paste0(args$sample_name, "\t", "NA"), fileConn)
    close(fileConn)
}
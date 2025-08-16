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
})

parser = argparse::ArgumentParser(description='Script to generate qc plots.')
parser$add_argument('cds_path', help='File with cds.')
parser$add_argument('umis_file', help='File with umis_per_barcode.')
parser$add_argument('sample_name', help='Sample name.')
parser$add_argument('empty_drops', help='RDS file from emptyDrops.')
parser$add_argument('--specify_cutoff', type='integer', help='Optional. Specifies a cutoff rather than choosing a UMI cutoff automatically.')
args = parser$parse_args()

sample_name <- args$sample_name

# Read empty drops RDS object
emptydrops_data <- readRDS(args$empty_drops)


#################################################
##   Functions for manuipulating cds objects   ##
##          performing well checks             ##
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

    # Get ligation plates
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

        # Map each ligation barcode to its plate
    cds$Ligation_plate <- vapply(lig_nums, function(n) {
      name <- names(lig_info)[sapply(lig_info, function(x) n %in% x)]
      if (length(name) == 0) return(NA_character_) else return(name)
    }, character(1))

    return (cds)
}

# Function to return labels/text for stat_summary
# Returns labels for ymax data (here, it is max number of cells)
# n_fun <- function(y){return(data.frame(y=max(y), label = paste0(length(y))))}


# Extract RT barcode and plate information from cell name
# Write csv for UMI and mitochondrial UMI summary stats by RT barcode
# Create mito and UMI plot for each RT barcode
# returns a list with a cds with the extracted information, mito and umi plots 
rt_stats <- function(sample_name, cds) {
  temp_cds <- extractBarcode(cds)
   
  mito_rt_plot <- ggplot(data.frame(colData(temp_cds)), aes(x=RT_barcode, y=perc_mitochondrial_umis)) +
    geom_violin(aes(fill="aquamarine")) +
    geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
    theme_light() +
    theme(axis.text.x=element_blank(),
          text=element_text(size=14)) +
    geom_hline(yintercept = 10, linetype="dotted", ) +
    scale_y_continuous(limits=c(0,max(cds$perc_mitochondrial_umis) + 5)) +
    xlab("") +
    ylab("% Mito UMIs") +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, vjust = -0.3) +
    theme(legend.position="none")


  # Num of UMIs by RT barcode 
  umi_rt_plot <- ggplot(data.frame(colData(temp_cds)), aes(x=RT_barcode, y=n.umi)) +
  # facet_wrap(~sample, nrow=1, drop=FALSE, scales="free_x") +
    geom_violin(aes(fill="salmon")) +
    geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
    # theme_bbi() +
    theme_light() +
    theme(axis.text.x=element_blank(),
        text=element_text(size=14)) +
    scale_y_log10() + 
    xlab("") +
    ylab("UMIs") +
    theme(legend.position="none")
  
  data_list <- list(temp_cds, mito_rt_plot, umi_rt_plot)
  return (data_list)
    
}


##### UMI Violin plots ######

# Generate UMIs plots by RT plate, ligation plate and p5 plate
# Return a combined UMI plot
plot_umis <- function(cds) {
  n_fun <- function(y){return(data.frame(y=max(y), label = paste0(length(y))))}
  print("here 121")

  cds <- extractBarcode(cds)
  print(cds)

  umi_plot_rt <- ggplot(data.frame(colData(cds)), aes(x=RT_barcode, y=n.umi)) +
  # facet_wrap(~plate, nrow=1, drop=TRUE, scales="free_x") +
  geom_violin(fill="salmon") +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light(base_size=10) +
  scale_y_log10(labels=trans_format("log10", math_format(10^.x)),
                limits=c(100, max(cds$n.umi) *30 )) +
  xlab("RT_plate") +
  ylab("UMIs") +
  stat_summary(fun.data = n_fun, geom = "text", angle=90, hjust = -0.1, size=3) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y= element_text(margin=margin(r=10))) 

  # barcode_list <- list("RT_plate", "Ligation_plate", "P5_barcode")
  # umi_plot_list <- list()

  # plot_list <- lapply(barcode_list, function(barcode) {

  #   barcode_name <- gsub("_", " ", barcode)

  #   print("here 130")
  #   print(head(data.frame(colData(cds))))
  #   print("barcode:")
  #   print(barcode)
  #   # Plot UMIs by RT, Ligation and P5 plates 
  #   umi_plot <- ggplot(data.frame(colData(cds)), aes(x=barcode, y=n.umi)) +
  #     geom_violin(fill="salmon")+
  #     geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  #     theme_light(base_size=10) +
  #     scale_y_log10(labels=trans_format("log10", math_format(10^.x)),
  #                   limits=c(100, max(cds$n.umi) *30 )) +
  #     xlab(barcode_name) +
  #     ylab("UMIs") +
  #     stat_summary(fun.data = n_fun, geom = "text", angle=90, hjust = -0.1, size=3) +
  #     theme(legend.position="none") +
  #     theme(axis.text.x = element_text(angle=90, hjust=1),
  #           axis.title.x = element_text(margin=margin(t=10)),
  #           axis.title.y= element_text(margin=margin(r=10))) 

  #     print("here 149")

  #     # umi_plot_list[[barcode]] <- umi_plot
  #   })

  # print("here 152")
  # flat_plot_list <- unlist(umi_plot_list, recursive=FALSE)

  # combined_umi <- plot_grid(flat_plot_list$RT_plate,
  #                           plot_spacer(),
  #                           flat_plot_list$Ligation_plate,
  #                           plot_spacer(), 
  #                           flat_plot_list$P5_barcode,
  #                           align="hv", 
  #                           ncol=1,
  #                           rel_heights = c(1, 0.25, 1, 0.25, 1))

  print("here 159")

  return (combined_umi)
}


#################################################
##       Function for Pseudobulk Plots         ##
#################################################

### Get average gene expression for each sample and RT barcode 
calc_pseudobulk <- function(cds) {
  
  barcode_list <- list("RT_barcode", "Ligation_plate", "P5_barcode")
  
  norm_counts <- normalized_counts(cds_filt)   
  plot_list = list()
  
  plot_list <- lapply(barcode_list, function(barcode) {
    print(barcode)
    
    psb_list <- list()
    
    # Precompute barcode groupings
    barcode_vec <- as.character(cds_filt[[barcode]])
    head(cds_filt[[barcode]])
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
    rownames(pseudo_bulkDF) <- paste(rowData(cds_filt)$gene_short_name, rowData(cds_filt)$id, sep="_") ## Set gene names as rownames
    
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
      # scale_fill_gradient(low = "darkblue",high = "yellow",
      #                     limits = c(0, 1.000)) + 
      theme_minimal(base_size=8) +
      xlab(colnames(psb)[1]) + 
      ylab(colnames(psb)[2]) +
      heatmap_theme 
    
    heatmap
    # scale_fill_continuous(low = "yellow", high = "purple", limits = c(0.5, 1.000))'
    
    plot_list[[paste0("heatmap_",barcode)]] <- heatmap
    
    # ggsave(paste0(barcode, "_pseudobulk_heatmap.png"), heatmap)
    
    
    ### Create distribution of barcode correlations
    hist_plot <- ggplot(psb, aes(x = .data[["pearson correlation"]])) +
      geom_histogram(binwidth = .003, fill = "grey", color = "black") +
      theme_bw(base_size=12) +
      # coord_cartesian(xlim = c(min(psb$value), 1)) +
      coord_cartesian(xlim = c(.60, 1)) +
      geom_vline(xintercept = 0.90, color = "red", linetype = "dashed", linewidth = 0.5) + 
      ggtitle(barcode_name)
    
    hist_plot
    
    plot_list[[paste0("hist_",barcode)]] <- hist_plot
    
    # ggsave(paste0(barcode, "_pseudobulk_hist_plot.png"), hist_plot)
    
    return(plot_list)
    
    # write.table(psb, paste0(barcode, "_psuedobulk_correlations.txt"), sep="\t", quote=FALSE, row.names = FALSE)
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
  

  # ggsave("combined_heatmaps.png", combined_heatmap, width=5, height=10, dpi=800, units='in')
  
  
  hist_plots <- unlist(all_plots[2], recursive=FALSE)
  hist_plots
  
  combined_hist <- plot_grid(flat_plot_list$hist_RT_barcode, 
                                 plot_spacer(),
                                 flat_plot_list$hist_Ligation_plate,
                                 plot_spacer(), 
                                 flat_plot_list$hist_P5_barcode,
                                 ncol=1, 
                                 rel_heights = c(1, 0.25, 1, 0.25, 1),
                                 align="hv")
  
  combined_hist
  
  # ggsave("combined_histograms.png", combined_hist, width=5, height=10, dpi=800, units='in')
  
}






#################################################
##         Function to generate UMAPs          ##
#################################################

plot_cells_simp <- function(cds,
                            colpal,
                            x=1,
                            y=2,
                            reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                            color_cells_by="cluster",
                            group_cells_by=c("cluster", "partition"),
                            trajectory_graph_color="grey28",
                            trajectory_graph_segment_size=0.75,
                            norm_method = c("log", "size_only"),
                            label_cell_groups = TRUE,
                            label_groups_by_cluster=TRUE,
                            group_label_size=2,
                            labels_per_group=1,
                            graph_label_size=2,
                            cell_size=0.35,
                            cell_stroke= I(cell_size / 2),
                            alpha = 1,
                            min_expr=0.1,
                            rasterize=FALSE) {
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimensions with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >=max(x,y),
                          msg = paste("x and/or y is too large. x and y must",
                                      "be dimensions in reduced dimension",
                                      "space."))
  if(!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                  "pseudotime") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must one of",
                                        "'cluster', 'partition', 'pseudotime,",
                                        "or a column in the colData table."))
    
    if(color_cells_by == "pseudotime") {
      tryCatch({pseudotime(cds, reduction_method = reduction_method)},
               error = function(x) {
                 stop(paste("No pseudotime for", reduction_method,
                            "calculated. Please run order_cells with",
                            "reduction_method =", reduction_method,
                            "before attempting to color by pseudotime."))})
      
    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))
  
  norm_method = match.arg(norm_method)
  group_cells_by=match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  
  if (rasterize){
    plotting_func <- ggrastr::geom_point_rast
  }else{
    plotting_func <- ggplot2::geom_point
  }
  
  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[,c(x,y)])
  
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster"){
    data_df$cell_group <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    data_df$cell_group <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else{
    stop("Error: unrecognized way of grouping cells.")
  }
  
  if (color_cells_by == "cluster"){
    data_df$cell_color <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "pseudotime") {
    data_df$cell_color <-
      tryCatch({pseudotime(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]}, error = function(e) {NULL})
  } else{
    data_df$cell_color <- colData(cds)[data_df$sample_name,color_cells_by]
  }
  
  
  
  if (label_cell_groups && is.null(color_cells_by) == FALSE){
    if (is.null(data_df$cell_color)){
      if (is.null(genes)){
        message(paste(color_cells_by, "not found in colData(cds), cells will",
                      "not be colored"))
      }
      text_df = NULL
      label_cell_groups = FALSE
    }else{
      if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
        
        if (label_groups_by_cluster && is.null(data_df$cell_group) == FALSE){
          text_df = data_df %>%
            dplyr::group_by(cell_group) %>%
            dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
            dplyr::group_by(cell_color, .add=TRUE) %>%
            dplyr::mutate(per=dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_group) %>%
            dplyr::top_n(labels_per_group, per)
        } else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>%
            dplyr::mutate(per=1)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_color) %>%
            dplyr::top_n(labels_per_group, per)
        }
        
        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
      } else {
        message(paste("Cells aren't colored in a way that allows them to",
                      "be grouped."))
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }
  
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
    
    if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        g <- g + geom_point(color=I("gray"), size=I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE,
                            alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't",
                      "color cells by cluster"))
      } else{
        g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE,
                            alpha = alpha)
      }
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    } else if (class(data_df$cell_color) == "numeric"){
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by, option="C")
    } else {
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }
    

  
#  if(label_cell_groups) {
#    g <- g + ggrepel::geom_text_repel(data = text_df,
#                                      mapping = aes_string(x = "text_x",
#                                                           y = "text_y",
#                                                           label = "label"),
#                                      size=I(group_label_size))
#      g <- g + theme(legend.position="none")
#  }
  
  g <- g +
    monocle3:::monocle_theme_opts() +
    xlab(paste(reduction_method, x)) +
    ylab(paste(reduction_method, y)) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) + 
    scale_color_manual(values=colpal) 
  g
}

#################################################
##             Start QC checks                 ##
#################################################


gen_plots <- function(sample_name, sample_path) {
  
  # samp_cds <- NULL
  # if (file.size(paste0(sample_path,"/bpcells_matrix_dir/col_names")) != 0) {
  #   samp_cds <- load_monocle_objects(sample_path)
  # } else {
  #   samp_cds <- readRDS(paste0(sample_path,"/cds_object.rds"))
  # }

  # garnett_mods <- names(colData(samp_cds))[grepl("garnett_type", names(colData(samp_cds)))]  
 

  samp_cds <- tryCatch({

    samp_cds <- load_monocle_objects(sample_path)
    print(samp_cds)

    print(head(colData(samp_cds)))
    # Generate UMI and mitochondrial stats by rt barcode 
    cds_rt_data <- rt_stats(sample_name, samp_cds) # Returns a list of data with cds and list of plots 
    print("here at 593")
    samp_cds <- cds_rt_data[[1]]

    # garnett_mods <- names(colData(samp_cds))[grepl("garnett_type", names(colData(samp_cds)))]  
    # samp_cds <- extractBarcode(samp_cds)

    # Generate UMI plots by RT, Ligation and P5 plates
    print("here 641")
    umi_plot <- plot_umis(samp_cds)

    print("here 643")
    # If Empty drops was ran, then filter out cds object for Empty Drops FDR < 0.01
    if(is(emptydrops_data, 'DFrame')) {
      keep_na <- samp_cds[,is.na(colData(samp_cds)$emptyDrops_FDR)]
      keep_cells <- samp_cds[,!is.na(colData(samp_cds)$emptyDrops_FDR) & colData(samp_cds)$emptyDrops_FDR <= 0.01]
      samp_cds <- combine_cds(list(keep_na, keep_cells), sample_col_name="og_cds")
    }

    # samp_cds <- detect_genes(samp_cds)

    # calc_pseudobulk(cds_filt)

    ggsave(paste0(sample_name, "_umi.png"), umi_plot, width=5, height=8, dpi=600, units="in")
    # ggsave(paste0(sample_name, "_pseudobulk_heatmap.png"), PLOT, width=5, height=10, dpi=800, units='in')
    # ggsave(paste0(sample_name, "_pseudobulk_histogram.png"), PLOT, width=5, height=10, dpi=800, units='in')
    
  }, error = function(e) {

     file_name <- paste0(sample_name, "_umi.png")
     ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for UMI plot")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="UMI", y = "barcode")
     ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')

    file_name <- paste0(sample_name, "_pseudobulk_heatmap.png")
    ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for Pseudobulk correlations")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none")
    ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')

    file_name <- paste0(sample_name, "_pseudobulk_histogram.png")
    ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for Pseudobulk correlations")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none")
    ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')

  })

  plot <- ggplot(as.data.frame(pData(samp_cds)), aes(n.umi, perc_mitochondrial_umis)) +
    geom_point(size=.5) + theme_bw() +
    lims(x=c(0, 100)) +
    labs(x="# of UMIs", y="Percent of UMIs from mitochondrial genome") +
    geom_rug(size=0.1) +
    scale_x_log10() + geom_density(aes(n.umi, (..density.. * max(pData(samp_cds)$perc_mitochondrial_umis))))
    ggsave(paste0( sample_name, '_cell_qc.png'), plot = plot, units = "in", width = 3.5*1.3, height = 3.5)
  
    samp_cds
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
      geom_hline(aes(yintercept = cutoff), linetype="dotted", size = .5, color = "firebrick2") + 
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
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
    cds$plate <- sapply(strsplit(as.character(meta$RT_barcode), "-"), `[`, 1)
    return (cds)
}

# Function to return labels/text for stat_summary
# Returns labels for ymax data (here, it is max number of cells)
n_fun <- function(y){return(data.frame(y=max(y), label = paste0(length(y))))}


# Extract RT barcode and plate information from cell name
# Write csv for UMI and mitochondrial UMI summary stats by RT barcode
# Create mito and UMI plot for each RT barcode
# returns a list with a cds with the extracted information, mito and umi plots 
rt_stats <- function(sample_name, cds) {
  temp_cds <- extractBarcode(cds)


  umi_rt_stats <- data.frame(colData(temp_cds)) %>% 
                group_by(RT_barcode) %>%
                dplyr::summarise(min=min(n.umi),
                                 q1=quantile(n.umi, probs=c(0.25)),
                                 med=median(n.umi),
                                 mean=mean(n.umi),
                                 q3=quantile(n.umi, probs=c(0.75)),
                                 max=max(n.umi)) %>%
                data.frame()

  write.csv(umi_rt_stats, file=paste0(sample_name, "_umi_rt_stats.csv"), quote=FALSE, row.names=FALSE)


  mito_rt_stats <- data.frame(colData(temp_cds)) %>% 
    group_by(RT_barcode) %>%
    summarise(min=min(perc_mitochondrial_umis),
              q1=quantile(perc_mitochondrial_umis, probs=c(0.25)),
              med=median(perc_mitochondrial_umis),
              mean=mean(perc_mitochondrial_umis),
              q3=quantile(perc_mitochondrial_umis, probs=c(0.75)),
              max=max(perc_mitochondrial_umis)) %>%
    data.frame()

  write.csv(mito_rt_stats, file=paste0(sample_name, "_mito_rt_stats.csv"), quote=FALSE, row.names=FALSE)
   
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


# 
recovery_plots <- function(cds) {
    # Violin plot of percent mitochondrial umis by RT barcode with threshold line at 10% 


}


# Perform well check for each RT well 
# Looks at number of UMIs, percent mitochondrial UMIs and 
# percent of a UMAP cluster in a well 
well_check <- function(sample_name, cds) {

  # Add UMAP coordinates to the colData for easy plotting outside of monocle3 
  cds$UMAP1 <- reducedDim(cds, "UMAP")[,1]
  cds$UMAP2 <- reducedDim(cds, "UMAP")[,2]

  # Extract the meta info into a data frame 
  meta <- data.frame(colData(cds))
 
  # Add cluster information to colData for each cell 
  colData(cds)$clusters = clusters(cds)

  # Generate random colors for easy viewing in well checks and umaps 
  nColor <- length(levels(cds$clusters))
  colpal <- randomcoloR::distinctColorPalette(k = nColor)
  # pie(rep(1, nColor), col = colpal)


  # Calculate percentages of each cluster 
  meta <- data.frame(colData(cds))
  clusterCounts <- meta %>% 
    dplyr::group_by(RT_barcode) %>%
    dplyr::count(clusters) %>%
    data.frame()

 # Number of cells for each RT barcode and cluster
  rt_key <- data.frame(table(meta$RT_barcode))
  cluster_key <- data.frame(table(meta$clusters))

  # Total number of cells for each RT barcode 
  clusterCounts$RTtot <- rt_key[match(clusterCounts$RT_barcode, rt_key$Var1), 'Freq']

  # Total number of cells per cluster / total RT barcode cells
  clusterCounts$perRT <- clusterCounts$n/clusterCounts$RTtot * 100

  # Total number of cells per cluster / total number of cells in sample
  clusterCounts$perSamp <- clusterCounts$n/nrow(colData(cds)) * 100

  # Find overall percent of each cell in each cluster by RT barcode 
  clusterCounts$clustTot <- cluster_key[match(clusterCounts$clusters, cluster_key$Var1), 'Freq']
  clusterCounts$perClust <- clusterCounts$n/clusterCounts$clustTot * 100

  clusterCounts$plate <- substr(clusterCounts$RT_barcode, 1, 3)
  meta$plate <- substr(meta$RT_barcode, 1, 3)


  # UMAP of all data to use as background cluster
  bg_data <- data.frame(colData(cds))[,c('UMAP1', 'UMAP2')]

  # # Violin plot of percent mitochondrial umis by RT barcode with threshold line at 10% 
  # mito_rt_plot <- ggplot(meta, aes(x=RT_barcode, y=perc_mitochondrial_umis)) +
  #   geom_violin(aes(fill="aquamarine")) +
  #   geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  #   theme_light() +
  #   theme(axis.text.x=element_blank(),
  #         text=element_text(size=14)) +
  #   geom_hline(yintercept = 10, linetype="dotted", ) +
  #   scale_y_continuous(limits=c(0,max(cds$perc_mitochondrial_umis) + 5)) +
  #   xlab("") +
  #   ylab("% Mito UMIs") +
  #   stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, vjust = -0.3) +
  #   theme(legend.position="none")


  # # Num of UMIs by RT barcode 
  # umi_rt_plot <- ggplot(meta, aes(x=RT_barcode, y=n.umi)) +
  # # facet_wrap(~sample, nrow=1, drop=FALSE, scales="free_x") +
  #   geom_violin(aes(fill="salmon")) +
  #   geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  #   # theme_bbi() +
  #   theme_light() +
  #   theme(axis.text.x=element_blank(),
  #       text=element_text(size=14)) +
  #   scale_y_log10() + 
  #   xlab("") +
  #   ylab("UMIs") +
  #   theme(legend.position="none")

  # Percent cluster of total sample 
  perc_clust_plot <- ggplot(clusterCounts, aes(x=RT_barcode, y=perSamp, fill=clusters)) +
    geom_bar(stat="identity") +
    theme_light() +
    scale_fill_manual(values=colpal) +
    theme(axis.text.x=element_blank(),
          text=element_text(size=14)) +
    xlab("") +
    ylab("% Total Sample") +
    theme(legend.position="none")


  # Percent cluster of RT well 
  perc_clust_rt <- ggplot(clusterCounts, aes(x=RT_barcode, y=perRT, fill=clusters)) +
    geom_bar(stat="identity") +
    theme_light() +
    scale_fill_manual(values=colpal) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          text=element_text(size=14)) +
    xlab("RT Wells") +
    ylab("% of RT Well") +
    theme(legend.position="bottom")

  # all 3 well check plots combined 
  # well_check_combined <- plot_grid(mito_rt_plot, umi_rt_plot, 
  #                        perc_clust_plot, perc_clust_rt, 
  #                        nrow=4, align='hv',
  #                        rel_heights = c(1,1,1,2))


  # ggsave(paste0(sample_name, "_wellcheck.png"), well_check_combined, width=15, height = 30)

  plot_list <- list(perc_clust_plot, perc_clust_rt, colpal)

  return (plot_list) # return colors used for clusters and list of well check plots
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
  
  samp_cds <- readRDS(sample_path)

  garnett_mods <- names(colData(samp_cds))[grepl("garnett_type", names(colData(samp_cds)))]  

  samp_cds <- tryCatch({

    # If Empty drops was ran, then filter out cds object for Empty Drops FDR < 0.01
    if(is(emptydrops_data, 'DFrame')) {
      keep_na <- samp_cds[,is.na(colData(samp_cds)$emptyDrops_FDR)]
      keep_cells <- samp_cds[,!is.na(colData(samp_cds)$emptyDrops_FDR) & colData(samp_cds)$emptyDrops_FDR <= 0.01]
      samp_cds <- combine_cds(list(keep_na, keep_cells), sample_col_name="og_cds")
    }

    # Generate UMI and mitochondrial stats by rt barcode 
    cds_rt_data <- rt_stats(sample_name, samp_cds) # Returns a list of data with cds and list of plots 
    samp_cds <- cds_rt_data[[1]]
    # ggsave("umi_plot", rt_plots[[3]])


    # Subsample cells from cds object if > 300,000 cells for preprocessing and clustering
    # Original cds will be used for RT umi and mitochondrial stats

    sub_cds <- NULL 
    wellcheck_data <- NULL 

    if(nrow(colData(samp_cds)) >300000) {

      sub_cds <- samp_cds[, samp_cds$cell %in% sample(samp_cds$cell, 300000)]
      sub_cds <- preprocess_cds(sub_cds)
      sub_cds <- reduce_dimension(sub_cds)
      sub_cds <- cluster_cells(sub_cds, k=ceiling(sqrt(dim(sub_cds)[2])*0.25))

      # Reduce number of clusters if clusters are > 12 
      if (dim(table(clusters(sub_cds))) > 12 ) {
        sub_cds <- cluster_cells(sub_cds, k=ceiling(sqrt(dim(sub_cds)[2])*0.75))
      }  

      # Increase number of clusters if clusters are < 5
      if (dim(table(clusters(sub_cds))) < 5 ) {
        sub_cds <- cluster_cells(sub_cds, k=ceiling(sqrt(dim(sub_cds)[2])*0.1))
      }  

      wellcheck_data = well_check(sample_name, sub_cds)
      colpal = wellcheck_data[[3]]

      # Plot umap 
      file_name <- paste0(sample_name, "_UMAP.png")
      ggp_obj <- suppressMessages(plot_cells_simp(sub_cds, colpal) + theme(text = element_text(size = 8)))
      ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')

    } else {

      samp_cds <- preprocess_cds(samp_cds)
      samp_cds <- reduce_dimension(samp_cds)
      samp_cds <- cluster_cells(samp_cds, k=ceiling(sqrt(dim(samp_cds)[2])*0.25))

      # Reduce number of clusters if clusters are > 12 
      if (dim(table(clusters(samp_cds))) > 12 ) {
        samp_cds <- cluster_cells(samp_cds, k=ceiling(sqrt(dim(samp_cds)[2])*0.75))
      }  

      # Increase number of clusters if clusters are < 5
      if (dim(table(clusters(samp_cds))) < 5 ) {
        samp_cds <- cluster_cells(samp_cds, k=ceiling(sqrt(dim(sub_cds)[2])*0.1))
      }  

      wellcheck_data = well_check(sample_name, samp_cds)
      colpal = wellcheck_data[[3]]
      
#    file_name <- paste0(sample_name, "_UMAP.png")
#    png(file_name, width = 5, height = 5, res = 600, units = "in")
#    print(suppressMessages(plot_cells_simp(samp_cds, colpal) + theme(text = element_text(size = 8))))
#    dev.off()

      # Plot umap
      file_name <- paste0(sample_name, "_UMAP.png")
      
      ggp_obj <- suppressMessages(plot_cells_simp(samp_cds, colpal) + theme(text = element_text(size = 8)))
      ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')
      
    }

    # Combine mito, umi, and rt well check plots 
    well_check_combined <- plot_grid(cds_rt_data[[2]], cds_rt_data[[3]], # mito and umi rt barcode plots 
                    wellcheck_data[[1]], wellcheck_data[[2]], # perc cluster plots from rt well check
                    nrow=4, align='hv',
                    rel_heights = c(1,1,1,2))


    ggsave(paste0(sample_name, "_wellcheck.png"), well_check_combined, width=15, height = 30)


    for (mod in garnett_mods) {
#      file_name <- paste0(sample_name, "_", gsub("garnett_type_", "", mod) ,"_Garnett.png")
#      png(file_name, width = 7, height = 5, res = 600, units = "in")
#      print(suppressMessages(plot_cells_simp(samp_cds, color_cells_by = mod) + theme(text = element_text(size = 8)) + theme(legend.position = "none")))
#      dev.off()

      file_name <- paste0(sample_name, "_", gsub("garnett_type_", "", mod) ,"_Garnett.png")
      ggp_obj <- suppressMessages(plot_cells_simp(samp_cds, color_cells_by = mod) + theme(text = element_text(size = 8)) + theme(legend.position = "none"))
      ggsave(filename=file_name, ggp_obj, device='png', width=7, height=5, dpi=600, units='in')
    }
    
    samp_cds
  }, error = function(e) {

#    png(paste0(sample_name, "_UMAP.png"), width = 5, height = 5, res = 600, units = "in")
#    print(ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for UMAP")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="UMAP 1", y = "UMAP 2"))
#    dev.off()

     file_name <- paste0(sample_name, "_UMAP.png")
     ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for UMAP")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="UMAP 1", y = "UMAP 2")
     ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')
   
#    png(paste0(sample_name, "_wellcheck.png"), width = 5, height = 5, res = 600, units = "in")
#    print(ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for WellCheck")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="RT Barcodes", y = "UMIs"))
#    dev.off()

    file_name <- paste0(sample_name, "_wellcheck.png")
    ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for WellCheck")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="RT Barcodes", y = "UMIs")
     ggsave(filename=file_name, ggp_obj, device='png', width=5, height=5, dpi=600, units='in')

    write.table("Insufficent UMIs for RT stats", file=paste0(sample_name, "_mito_rt_stats.csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table("Insufficent UMIs for RT stats", file=paste0(sample_name, "_umi_rt_stats.csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)

    for (mod in garnett_mods) {
#      png(paste0(sample_name, "_", gsub("garnett_type_", "", mod) ,"_Garnett.png"), width = 7, height = 5, res = 600, units = "in")
#      print(ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for UMAP")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="UMAP 1", y = "UMAP 2"))
#      dev.off()

      file_name <- paste0(sample_name, "_", gsub("garnett_type_", "", mod) ,"_Garnett.png")
      ggp_obj <- ggplot() + geom_text(aes(x = 1, y = 1, label = "Insufficient cells for UMAP")) + monocle3:::monocle_theme_opts() + theme(legend.position = "none") + labs(x="UMAP 1", y = "UMAP 2")
      ggsave(filename=file_name, ggp_obj, device='png', width=7, height=5, dpi=600, units='in')

    }
    samp_cds  
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
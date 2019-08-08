#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(argparse)
  library(shiny)
  library(stringr)
  library(monocle3)
})

parser = argparse::ArgumentParser(description='Script to generate experiment dashboard.')
parser$add_argument('input_folder', help='input folder.')
parser$add_argument('--all_dups', required=TRUE, help='all dup file')
args = parser$parse_args()

output_folder <- args$input_folder
sample_folds <- list.files(output_folder)
dedup <- readLines(args$all_dups)
dedup <- dedup[seq(2, length(dedup), by=2)]
dedup_df <- as.data.frame(stringr::str_split_fixed(dedup, " +", 4))
dedup_df$sample <- stringr::str_split_fixed(dedup_df$V1, ":", 2)[,2]

project_name <- unlist(stringr::str_split(output_folder, "/"))
project_name <- project_name[[length(project_name)]]

barn <- ""
sent <- ""

top <- HTML('<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
    <style>
        .navbar-brand
        {
          font-size: xx-large !important;
          display: flex;
          align-items: center;
          text-align: center;
        }
        .sidebar {
          position: fixed;
          top: 0;
          bottom: 0;
          left: 10px;
          z-index: 100; /* Behind the navbar */
          padding: 100px 0 0; /* Height of navbar */
          box-shadow: inset -1px 0 0 rgba(0, 0, 0, .1);
        }
        </style>
  </head>')


if ("Barnyard" %in% sample_folds) {
  cds <- readRDS(paste0(output_folder, "/Barnyard/Barnyard_cds.RDS"))
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
  
  ggsave("exp_dash/img/Barnyard_UMIs.png", plot = plot, units = "in", width = 3.5*1.3, height = 3.5)

  collision_rate <- round(sum(pData(cds)$collision/nrow(pData(cds))) * 200, 1)
  system(paste0("cp ", output_folder, "/Barnyard/knee_plot.png", " exp_dash/img/Barnyard_knee_plot.png"))
 barn <- list(HTML('<div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom">
                <h1 class="h3" id="rt">Barnyard</h1>
            </div>
          <nav>
              <div class="nav nav-tabs" id="navbarn-tab" role="tablist">
                <a class="nav-item nav-link active" id="navbarn-knee-tab" data-toggle="tab" href="#navbarn-knee" role="tab" aria-controls="navbarn-knee" aria-selected="true">Knee plot</a>
                <a class="nav-item nav-link" id="navbarn-barn-tab" data-toggle="tab" href="#navbarn-barn" role="tab" aria-controls="navbarn-barn" aria-selected="false">Barnyard</a>
                <a class="nav-item nav-link" id="navbarn-stats-tab" data-toggle="tab" href="#navbarn-stats" role="tab" aria-controls="navbarn-stats" aria-selected="false">Stats</a>
              </div>
          </nav>
          <div class="tab-content" id="nav-tabContent">
                <div class="tab-pane fade show active" id="navbarn-knee" role="tabpanel" aria-labelledby="navbarn-knee-tab">
                  <img src="img/Barnyard_knee_plot.png" width = 80%; class="rounded mx-auto d-block" alt="...">
                </div>
                <div class="tab-pane fade" id="navbarn-barn" role="tabpanel" aria-labelledby="navbarn-barn-tab">
                    <img src="img/Barnyard_UMIs.png" width = 80%; class="rounded mx-auto d-block" alt="...">
                </div>
                <div class="tab-pane fade" id="navbarn-stats" role="tabpanel" aria-labelledby="navbarn-stats-tab">
                    <table class="table table-hover">
                        <thead>
                          <tr>
                            <th scope="col"></th>
                            <th scope="col">Barnyard</th>
                          </tr>
                        </thead>
                        <tbody>
                          <tr>
                            <th scope="row">Total reads</th>
                            <td>'), dedup_df[dedup_df$sample == "Barnyard",]$V2, HTML('</td>
                          </tr>
                          <tr>
                            <th scope="row">Total UMIs</th>
                            <td>'), dedup_df[dedup_df$sample == "Barnyard",]$V3, HTML('</td>
                          </tr>
                          <tr>
                            <th scope="row">Duplication rate</th>
                            <td>'), dedup_df[dedup_df$sample == "Barnyard",]$V4, HTML('</td>
                          </tr>
                          <tr>
                              <th scope="row">Collision rate</th>
                              <td>'), paste0(collision_rate, "%"), HTML('</td>
                            </tr>
                        </tbody>
                      </table>
                </div>
                </div>'))
}

if ("Sentinel" %in% sample_folds) {
  system(paste0("cp ", output_folder, "/Sentinel/knee_plot.png", " exp_dash/img/Sentinel_knee_plot.png"))
  sent <- list(HTML(
    '             <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom">
                  <h1 class="h3" id="sent">Sentinel</h1>
              </div>
            <nav>
                <div class="nav nav-tabs" id="navsent-tab" role="tablist">
                  <a class="nav-item nav-link active" id="navsent-knee-tab" data-toggle="tab" href="#navsent-knee" role="tab" aria-controls="navsent-knee" aria-selected="true">Knee plot</a>
                  <a class="nav-item nav-link" id="navsent-umap-tab" data-toggle="tab" href="#navsent-umap" role="tab" aria-controls="navsent-umap" aria-selected="false">UMAP</a>
                  <a class="nav-item nav-link" id="navsent-stats-tab" data-toggle="tab" href="#navsent-stats" role="tab" aria-controls="navsent-stats" aria-selected="false">Stats</a>
                </div>
            </nav>
            <div class="tab-content" id="nav-tabContent">
                  <div class="tab-pane fade show active" id="navsent-knee" role="tabpanel" aria-labelledby="navsent-knee-tab">
                    <img src="img/Sentinel_knee_plot.png" width = 80%; class="rounded mx-auto d-block" alt="...">
                  </div>
                  <div class="tab-pane fade" id="navsent-umap" role="tabpanel" aria-labelledby="navsent-umap-tab">
                      <img src="img/Sentinel_UMAP.png" width = 80%; class="rounded mx-auto d-block" alt="...">
                  </div>
                  <div class="tab-pane fade" id="navsent-stats" role="tabpanel" aria-labelledby="navsent-stats-tab">
                      <table class="table table-hover">
                          <thead>
                            <tr>
                              <th scope="col"></th>
                              <th scope="col">Sentinel</th>
                            </tr>
                          </thead>
                                                  <tbody>
                          <tr>
                            <th scope="row">Total reads</th>
                            <td>'), dedup_df[dedup_df$sample == "Sentinel",]$V2, HTML('</td>
                          </tr>
                          <tr>
                            <th scope="row">Total UMIs</th>
                            <td>'), dedup_df[dedup_df$sample == "Sentinel",]$V3, HTML('</td>
                          </tr>
                          <tr>
                            <th scope="row">Duplication rate</th>
                            <td>'), dedup_df[dedup_df$sample == "Sentinel",]$V4, HTML('</td>
                          </tr>
                        </tbody>
                        </table>
                  </div>
                  </div>'
  ))
}
samp_html <- list()
good_samp_list <- list()
for(samp in sample_folds) {
  samp_np <- gsub("\\.", "", samp)
  if (samp %in% c("Sentinel", "Barnyard")) {
    next
  }
  if(!file.exists(paste0(output_folder, "/", samp, "/knee_plot.png"))) {
    next
  }
  good_samp_list <- append(good_samp_list, samp)
  system(paste0("cp ", output_folder, "/", samp, "/knee_plot.png", " exp_dash/img/", samp, "_knee_plot.png"))
  umis_per_cell <- read.table(paste0(output_folder, "/", samp, "/umis_per_cell_barcode.txt"))
  samp_html <- append(samp_html,
                  as.character(tags$div(class="tab-pane fade", id=paste0("nav", samp_np), role="tabpanel", `aria-labelledby`=paste0("nav", samp_np, "-tab"),
                   HTML(paste0('<br>
                          <table class="table table-hover">
                              <thead>
                                <tr>
                                  <th scope="col"></th>
                                  <th scope="col">Total reads</th>
                                  <th scope="col">Total UMIs</th>
                                  <th scope="col">Duplication rate</th>
                                  <th scope="col">Cells with >100 UMIs</th>
                                  <th scope="col">Cells with >1000 UMIs</th>
                                </tr>
                              </thead>
                              <tbody>
                                <tr>
                                  <th scope="row">', samp, '</th>
                                  <td>', dedup_df[dedup_df$sample == samp,]$V2 ,'</td>
                                  <td>', dedup_df[dedup_df$sample == samp,]$V3 ,'</td>
                                  <td>', dedup_df[dedup_df$sample == samp,]$V4 ,'</td>
                                  <td>', sum(umis_per_cell$V3 > 100) ,'</td>
                                  <td>', sum(umis_per_cell$V3 > 1000) ,'</td>
                                </tr>
                              </tbody>
                            </table>
                    <img src="img/', samp, '_knee_plot.png" width = 80%; class="rounded mx-auto d-block" alt="...">'
                  
                  ))),
                  '</div>'))
}

good_samp_tabs <- list()
for(samp in good_samp_list) {
  samp_np <- gsub("\\.", "", samp)
  good_samp_tabs <- append(good_samp_tabs, 
                           as.character(tags$a(class="nav-item nav-link", 
                                  id=paste0("nav", samp_np, "-tab"), 
                                  href=paste0("#nav", samp_np), 
                                  `data-toggle`="tab",
                                  role="tab", 
                                  `aria-controls`=paste0("nav", samp_np), 
                                  `aria-selected`="false", 
                                  samp)))
}

body <- tags$body(
  list(
  HTML('  <nav class="navbar navbar-expand-md sticky-top navbar-light" style="background-color: #e3f2fd;">
        <div class="navbar-collapse collapse w-100 order-1 order-md-0 dual-collapse2">
            <ul class="navbar-nav mr-auto">
                <img src="img/BBI_Logo_Horizontal_RGB_Pantone Equivilant.png" height="70" class="d-inline-block align-top" alt="">
            </ul>
        </div>
        <div class="mx-auto order-0">
            <a class="navbar-brand mx-auto" href="#">'), paste('Experiment', project_name, 'QC Dashboard'), HTML('</a>
        </div>
        <div class="navbar-collapse collapse w-100 order-3 dual-collapse2">
        </div>
    </nav>
    <div class="container-fluid">
      <div class="row">
        <nav class="col-md-2 d-none d-md-block bg-light sidebar">
          <div class="sidebar-sticky">
            <ul class="nav flex-column">
              <li class="nav-item">
                <a class="nav-link active" href="#summary">
                  <span data-feather="home"></span>
                  Summary Statistics <span class="sr-only">(current)</span>
                </a>
              </li>
              <li class="nav-item">
                <a class="nav-link" href="#barn">
                  <span data-feather="file"></span>
                  Barnyard
                </a>
              </li>
              <li class="nav-item">
                <a class="nav-link" href="#sent">
                  <span data-feather="shopping-cart"></span>
                  Sentinel
                </a>
              </li>
              <li class="nav-item">
                <a class="nav-link" href="#sent">
                  <span data-feather="shopping-cart"></span>
                  Samples
                </a>
              </li>
            </ul>
          </div>
        </nav>
       <main role="main" class="col-md-9 ml-sm-auto col-lg-10 px-4" style="padding-top: 15px;">'),
  barn,
   sent,
  HTML('                <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom">
                  <h1 class="h3" id="sent">Samples</h1>
              </div>
        <nav>
                <div class="nav nav-tabs" id="navsamp1-tab" role="tablist">'
       ),
 HTML(unlist(good_samp_tabs)),
  HTML('</div>
            </nav>'),
 HTML('            <div class="tab-content" id="nav-tabContent">'),
  HTML(unlist(samp_html)),
    HTML('
    </main>
      </div>
      </div>
      
      
      <!-- Optional JavaScript -->
      <!-- jQuery first, then Popper.js, then Bootstrap JS -->
      <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1" crossorigin="anonymous"></script>
      <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js" integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM" crossorigin="anonymous"></script>
      ')
  )
)

end <- HTML('</html>')
  

fileConn<-file("exp_dash/exp_dash.html")
writeLines(c(as.character(top),as.character(body), as.character(end)), fileConn)
close(fileConn)





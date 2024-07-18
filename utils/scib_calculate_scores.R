suppressPackageStartupMessages({
  library(Seurat)
  library(funkyheatmap)
  library(tidyverse)
  library(scales)
  library(ggimage)
  library(cowplot)
  library(tibble)
  library(RColorBrewer)
  library(dynutils)
  library(stringr)
  library(Hmisc)
  library(plyr)
})

scib_calculate_scores <- function(scib_bm_table, Overall_Score_scale = F){
  bm <- scib_bm_table
  
  bio <- c('NMI_cluster.label', 'ARI_cluster.label', 'ASW_label', 'isolated_label_F1', 
           'isolated_label_silhouette', 'cLISI', 'cell_cycle_conservation')
  batch <- c('PCR_batch', 'ASW_label.batch', 'iLISI', 'graph_conn', 'kBET')
  # filter columns contain NA
  bio_keep <- c()
  for (x in bio) {
    if (!is.na(sum(bm[,x]))) {
      bio_keep <- c(bio_keep, x)
    }
  }
  batch_keep <- c()
  for (x in batch) {
    if (!is.na(sum(bm[,x]))) {
      batch_keep <- c(batch_keep, x)
    }
  }
  extract_tab <- bm[,c(bio_keep, batch_keep)]
  
  # calculate
  extract_tab <- extract_tab %>% 
    mutate(Batch_Correction = rowMeans(.[,batch_keep])) %>% 
    mutate(Bio_conservation = rowMeans(.[,bio_keep])) %>%
    mutate(Overall_Score = 0.4 * Batch_Correction + 0.6 * Bio_conservation) %>%
    dplyr::arrange(desc(Overall_Score))
  
  if (Overall_Score_scale) {
    extract_tab <- extract_tab %>% 
      dplyr::mutate(Overall_Score_scale = dynutils::scale_minmax(.$Overall_Score))
  }
  
  return(extract_tab)
}

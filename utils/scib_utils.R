
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

## SCIB data preprocessing

scib_calculate_scores <- function(scib_bm_table, Overall_Score_scale = F){
  bm <- scib_bm_table
  
  bio <- c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'isolated_label_F1', 
           'isolated_label_silhouette', 'cLISI', 'cell_cycle_conservation')
  batch <- c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn', 'kBET')
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


# data preprocessing from tab
scib_score2tab <- function(
    scib_summary_tab,
    Methods = "Methods",
    Features = "Features",
    Scaling = "Scaling",
    Output = "Output",
    Overall_Score = "Overall_Score",
    Bio_Score = "Bio_conservation",
    Batch_Score = "Batch_Correction"
){
  
  ### step 1: check columns names
  required_columns <- c(Methods, Features, Scaling, Output,
                        Overall_Score, Bio_Score, Batch_Score)
  missing_columns <- required_columns[!required_columns %in% colnames(scib_summary_tab)]
  if (length(missing_columns) > 0) {
    stop("The following required columns are missing in the input data: ", paste(missing_columns, collapse = ", "))
  }else{
    print("All required columns are present.")
  }
  
  data <- scib_summary_tab[required_columns]
  colnames(data) <- c("Methods", "Features", "Scaling", "Output",
                      "Overall_Score", "Bio_conservation", "Batch_Correction")
  # check output columns values
  check_column_value(data, column_name = "Output", c("Embedding", "Features", "Graph"))
  check_column_value(data, column_name = "Scaling", c("unscaled", "scaled"))
  check_column_value(data, column_name = "Features", c("HVG", "FULL"))
  
  # step 2: create ranking columns
  data <- data %>% 
    mutate(overall_rank = rank(-.$Overall_Score)) %>% 
    mutate(bio_rank = rank(-.$Bio_conservation)) %>%  
    mutate(batch_rank = rank(-.$Batch_Correction)) %>% 
    mutate_at(
      c("overall_rank", "bio_rank", "batch_rank"),
      function(x) {dynutils::scale_minmax(-x)}
    ) %>% 
    mutate(
      output_img = case_match(
        Output,
        "Features"  ~ "gene",
        "Embedding" ~ "embed",
        "Graph"     ~ "graph"
      )
    ) %>% 
    dplyr::select(
      Methods, Features, Scaling,
      overall_rank, bio_rank, batch_rank, 
      Output, output_img, 
      Overall_Score, Bio_conservation, Batch_Correction
    ) %>% 
    as.data.frame()
  
  print("Data preprocessing has been completed...")
  return(data)
}

# function for preparing scib_plotting data
scib_tab2plot <- function(
    preprocessed_tab, label_top3 = T, Image_path = "./images"
){
  
  # a list for collecting results
  plot_data = list()
  
  # preprocessed data from scib_preprocessing_bm
  data.plot <- preprocessed_tab
  
  # step 1: prepare row and columns info
  # row info
  row_info <- data.frame(id = data.plot$Methods)
  
  # create column_info
  column_info <- tibble(
    "id" = c("Methods", "output_img", "Features", "Scaling", "Overall_Score", "Bio_conservation", "Batch_Correction"),
    "id_colors" = c(NA, NA, NA, NA, "overall_rank", "bio_rank", "batch_rank"),
    "name" = c("Methods","Output",  "HVG", "Scaling", "Overall score", "Bioconserv.", "Batch correction"),
    "geom" = c("text", "image", "text", "text", "bar", "bar", "bar"),
    "group" = c("Text", "Text", "Image", "Text",
                "Overall Score", "Bio conservation", "Batch Correction"), # actually for palette
    "width" = c(8, 2, 2.5, 2, 2, 2, 2)
  )
  # set palettes for continuous variables
  # reference: https://matplotlib.org/stable/users/explain/colors/colormaps.html
  palettes <- list(
    "Bio conservation" = "YlGnBu",
    "Batch Correction" = "BuPu",
    "Overall Score" = "RdPu"
  )
  
  ### "id" includes elementals which are used to plot in the graph
  ### their relative values are in dataframe <scib_summary_plot>
  # methods : integration methods
  # features: HVG or FULL 
  # output_image: images for output types
  # scaling: data scaling or not
  # three scib summary indices
  
  ### "id_colors" includes colors applied to indices
  # NA means "black" (Default)
  # "bio_rank", "batch_rank", "overall_rank" are continuous variables from 0 to 1 (why to scale)
  
  ### "name" includes names for title in figure
  ### "geom" includes ggplot2 layer forms
  ### "group" for calculating position
  ### "width" indicates the custom widths used for drawing 
  # If width does not exist, it will be automatically completed.
  ### "overlay" for text overlay
  
  
  # step 2: initial plotting parameters
  # no point in making these into parameters
  row_height <- 1.1
  row_space <- .1
  row_bigspace <- .5
  col_width <- 1.1
  col_space <- .2
  col_bigspace <- .5
  
  # step 3: Determine row and column positions
  # row position
  if (!"group" %in% colnames(row_info) || all(is.na(row_info$group))) {
    row_info$group <- ""
    row_groups <- tibble(group = "")
    plot_row_annotation <- FALSE
  } else {
    plot_row_annotation <- TRUE
  }
  # define the background position
  row_pos <-
    row_info %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(group_i = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      row_i = row_number(),
      colour_background = group_i %% 2 == 1,
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      ysep = ifelse(do_spacing, row_height + 2 * row_space, row_space),
      y = - (row_i * row_height + cumsum(ysep)),
      ymin = y - row_height / 2,
      ymax = y + row_height / 2
    )
  
  # column position
  if (!"group" %in% colnames(column_info) || all(is.na(column_info$group))) {
    column_info$group <- ""
    plot_column_annotation <- FALSE
  } else {
    plot_column_annotation <- TRUE
  }
  
  # Custom widths are recommended
  column_info <-
    column_info %>%
    add_column_if_missing(width = col_width, overlay = FALSE) 
  
  # no need to change this part
  column_pos <-
    column_info %>%
    mutate(
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      xsep = case_when(
        overlay ~ c(0, -head(width, -1)),
        do_spacing ~ col_bigspace,
        TRUE ~ col_space
      ),
      xwidth = case_when(
        overlay & width < 0 ~ width - xsep,
        overlay ~ -xsep,
        TRUE ~ width
      ),
      xmax = cumsum(xwidth + xsep),
      xmin = xmax - xwidth,
      x = xmin + xwidth / 2
    )
  rownames(column_pos) = column_pos$id  
  # Even though it suggests that this feature has been removed, it actually still works!
  
  # step 4: create geom data
  # gather circle data
  ind_circle <- which(column_info$geom == "circle")
  if(length(ind_circle) > 0){
    dat_mat <- as.matrix(data[, ind_circle])
    col_palette <- data.frame(metric = colnames(dat_mat), 
                              group = column_info[match(colnames(dat_mat), column_info$id), "group"])
    
    col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
    
    circle_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                    function(x) rep(x, nrow(dat_mat)))), 
                              x0 = unlist(lapply(column_pos$x[ind_circle], 
                                                 function(x) rep(x, nrow(dat_mat)))), 
                              y0 = rep(row_pos$y, ncol(dat_mat)),
                              r = row_height/2*as.vector(sqrt(dat_mat))
    )
    for(l in unique(circle_data$label)){
      ind_l <- which(circle_data$label == l)
      circle_data[ind_l, "r"] <- rescale(circle_data[ind_l, "r"], to = c(0.05, 0.55), from = range(circle_data[ind_l, "r"], na.rm = T))
    }
    
    colors <- NULL
    
    
    for(i in 1:ncol(dat_mat)){
      palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
      colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
    }
    
    circle_data$colors <- colors
  }
  
  # gather bar data
  ind_bar <- which(column_info$geom == "bar")
  ind_bar <- column_info$id[ind_bar]  ## fixed bug
  dat_mat <- as.matrix(data.plot[, ind_bar])
  
  # define colors
  col_palette <- data.frame(metric = colnames(dat_mat), 
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])
  col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
  
  # define position 
  rect_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                function(x) rep(x, nrow(dat_mat)))),
                          method = rep(row_info$id, ncol(dat_mat)),
                          value = as.vector(dat_mat),
                          # fixed
                          xmin = unlist(lapply(column_pos[ind_bar, "xmin"] %>% as.vector %>% unlist(), 
                                               function(x) rep(x, nrow(dat_mat)))),
                          # fixed
                          xmax = unlist(lapply(column_pos[ind_bar, "xmax"]%>% as.vector %>% unlist(), 
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          xwidth = unlist(lapply(column_pos[ind_bar, "xwidth"], 
                                                 function(x) rep(x, nrow(dat_mat))))
  )
  rect_data <- rect_data %>%
    add_column_if_missing(hjust = 0) %>%
    mutate(
      xmin = xmin + (1 - value) * xwidth * hjust,
      xmax = xmax - (1 - value) * xwidth * (1 - hjust)
    )
  # you need to check the xmin and xmax
  
  # add color information
  colors <- NULL
  for(i in 1:ncol(dat_mat)){
    palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data.plot)-sum(is.na(dat_mat[,i])))
    colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
  }
  rect_data$colors <- colors
  
  
  # gather text data
  ind_text <- which(column_info$geom == "text")
  ind_text <- column_info$id[ind_text]
  dat_mat <- as.matrix(data.plot[, ind_text])
  
  text_data <- data.frame(label_value = as.vector(dat_mat), 
                          group = rep(colnames(dat_mat), each = nrow(dat_mat)),
                          xmin = unlist(lapply(column_pos[ind_text, "xmin"] %>% as.vector() %>% unlist(), 
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_text, "xmax"] %>% as.vector() %>% unlist(), 
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          size = 4, fontface = "plain", stringsAsFactors = F)
  
  text_data$colors <- "black"
  text_data[text_data$label_value == "HVG", "colors"] <- "darkgreen"
  text_data[text_data$label_value == "FULL", "colors"] <- "grey30"
  
  # replace scaled/unscaled with +/-
  text_data$label_value <- mapvalues(text_data$label_value, 
                                     from = c("scaled", "unscaled"), 
                                     to = c("+", "-"))
  
  text_data[text_data$label_value == "+" | text_data$label_value == "-", "size"] <- 5
  text_data[text_data$label_value == "+" | text_data$label_value == "-", "fontface"] <- "plain"
  
  # optional
  # ADD top3 ranking for each bar column
  if(label_top3){
    cols_bar <- unique(rect_data$label)
    cols_bar <- as.character(cols_bar[!is.na(cols_bar)])
    for(c in cols_bar){
      rect_tmp <- rect_data[rect_data$label == c,] %>% na.omit()
      rect_tmp <- add_column(rect_tmp, "label_value" = as.character(rank(-rect_tmp$value, ties.method = "min")))
      rect_tmp <- rect_tmp[rect_tmp$label_value %in% c("1", "2", "3"), c("label_value", "xmin", "xmax", "ymin", "ymax")]
      rect_tmp <- add_column(rect_tmp, "size" = 3, .after = "ymax")
      rect_tmp <- add_column(rect_tmp, "colors" = "black", .after = "size")
      rect_tmp <- add_column(rect_tmp, "fontface" = "plain", .after = "colors")
      rect_tmp <- add_column(rect_tmp, "group" = "top3", .after = "fontface")
      text_data <- bind_rows(text_data, rect_tmp)
    }
  }
  
  # add columns names
  df <- column_pos %>% filter(id != "Methods") %>% filter(id != "Ranking")
  colnames_ymax <- min(text_data$ymin) - 0.5
  colnames_ymin <- colnames_ymax - 3.75
  
  if (nrow(df) > 0) {
    text_data <-
      bind_rows(
        text_data,
        df %>% transmute(
          xmin = x, xmax = x, ymin = colnames_ymin, ymax = colnames_ymax,
          fontface = "plain", colors = "black",
          angle = 90, vjust = 0.5, hjust = 1,   # vjust = 1底部对齐, hjust = 1 右对齐；
          label_value = name, 
          size = 3.5
        )
      )
  }
  
  # generate row annotation
  if (plot_row_annotation) {
    row_annotation <-
      row_pos %>% 
      select(group, ymin, ymax) %>%
      group_by(group) %>%
      summarise(
        ymin = min(ymin),
        ymax = max(ymax),
        y = (ymin + ymax) / 2
      ) %>%
      ungroup() %>%
      mutate(xmin = -.5, xmax = 5) %>%
      filter(!is.na(group), group != "")
    
    text_data <- text_data %>% bind_rows(
      row_annotation %>%
        transmute(xmin, xmax, ymin = ymax + row_space, label_value = group %>% gsub("\n", " ", .), 
                  hjust = 0, vjust = .5, fontface = "bold", size = 4) %>%
        mutate(ymax = ymin + row_height)
    )
  }
  
  
  # gather image data
  ind_img_loc <- which(column_info$geom == "image")
  ind_img <- column_info$id[ind_img_loc]
  if(length(ind_img) > 0){
    dat_mat <- as.matrix(data.plot[, ind_img])
    
    images_paths <- c(
      paste0(Image_path, "/graph.png"),
      paste0(Image_path, "/embedding.png"),
      paste0(Image_path, "/matrix.png")
    )
    
    image_data <- data.frame(x = unlist(lapply(column_pos$x[ind_img_loc], 
                                               function(x) rep(x, nrow(dat_mat)))), 
                             y = rep(row_pos$y, ncol(dat_mat)),
                             image = mapvalues(dat_mat, from = c("graph", "embed", "gene"), 
                                               # change the dir path to your image path
                                               to = images_paths),
                             stringsAsFactors = FALSE
    )
    
  }
  
  suppressWarnings({
    minimum_x <- min(column_pos$xmin, text_data$xmin, na.rm = TRUE)
    maximum_x <- max(column_pos$xmax, text_data$xmax, na.rm = TRUE)
    minimum_y <- min(row_pos$ymin, text_data$ymin, na.rm = TRUE)
    maximum_y <- max(row_pos$ymax, text_data$ymax, na.rm = TRUE)
  })
  
  
  # step 5: create legend data 
  x_min_output <- minimum_x+0.5
  x_min_scaling <- minimum_x + 5.5
  x_min_ranking <- minimum_x + 10.5
  x_min_score <-  minimum_x + 17
  
  leg_max_y <- minimum_y - .5
  
  # Create legend for Output
  leg_min_x <- x_min_output
  output_title_data <- data.frame(xmin = leg_min_x, 
                                  xmax = leg_min_x+ 2, 
                                  ymin = leg_max_y - 1, 
                                  ymax = leg_max_y, 
                                  label_value = "Output", 
                                  hjust = 0, vjust = 0, 
                                  fontface = "bold",
                                  size = 3)
  
  output_img <- data.frame(x = leg_min_x+0.5,
                           y = c(leg_max_y-2, leg_max_y-3.2,leg_max_y-4.4),
                           image = images_paths
  )
  output_text <- data.frame(xmin = leg_min_x+1.5, 
                            xmax = leg_min_x+3, 
                            ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
                            ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
                            label_value = c("Gene", "Embed", "Graph"), 
                            hjust = 0, vjust = 0, 
                            fontface = "plain",
                            size = 3)
  text_data <- bind_rows(text_data, output_text, output_title_data)
  image_data <- bind_rows(image_data, output_img)
  
  
  # Create legend for scaling
  leg_min_x <- x_min_scaling
  scaling_title_data <- data.frame(xmin = leg_min_x, 
                                   xmax = leg_min_x+ 2, 
                                   ymin = leg_max_y - 1, 
                                   ymax = leg_max_y, 
                                   label_value = "Scaling", 
                                   hjust = 0, vjust = 0, 
                                   fontface = "bold",
                                   size = 3)
  
  scaling_text <- data.frame(xmin = c(leg_min_x, leg_min_x+1), 
                             xmax = c(leg_min_x+0.5, leg_min_x+3), 
                             ymin = c(rep(leg_max_y-2,2), rep(leg_max_y-3,2)), 
                             ymax = c(rep(leg_max_y-1,2), rep(leg_max_y-2,2)), 
                             label_value = c("+", " Scaled", "-", " Unscaled"), 
                             hjust = 0.5, vjust = 0, 
                             fontface = c("bold","plain", "bold", "plain"),
                             size = c(5,3,5,3))
  
  text_data <- bind_rows(text_data, scaling_title_data, scaling_text)
  
  # create legend for ranking colors
  leg_min_x <- x_min_ranking
  # rank_groups <- as.character(column_info[column_info$geom == "bar", "group"])
  rank_groups <- column_info[column_info$geom == "bar", "group"] %>% 
    as.vector() %>% 
    unlist() %>% 
    unname()
  rank_minimum_x <- list("Overall Score" = leg_min_x, 
                         "Batch Correction" = leg_min_x+1, 
                         "Bio conservation" = leg_min_x+2)
  leg_max_x <- leg_min_x+2
  rank_title_data <- data.frame(xmin = leg_min_x, 
                                xmax = leg_min_x+ 2, 
                                ymin = leg_max_y - 1, 
                                ymax = leg_max_y, 
                                label_value = "Ranking", 
                                hjust = 0, vjust = 0, 
                                fontface = "bold")
  
  for(rg in rank_groups){
    rank_palette <- colorRampPalette(rev(brewer.pal(9, palettes[[rg]])))(5)
    
    
    rank_data <- data.frame(xmin = rank_minimum_x[[rg]],
                            xmax = rank_minimum_x[[rg]] + .8,
                            ymin = seq(leg_max_y-4, leg_max_y - 2, by = .5),
                            ymax = seq(leg_max_y-3.5, leg_max_y -1.5, by = .5),
                            border = TRUE,
                            colors = rank_palette
    )
    rect_data <- bind_rows(rect_data, rank_data)
    
  }
  
  # create arrow for ranking
  arrow_data <- data.frame(x = leg_max_x + 1.5, 
                           xend = leg_max_x +1.5, 
                           y = leg_max_y-4, 
                           yend = leg_max_y -1.5)
  
  # add text next to the arrow
  arrow_text <- data.frame(xmin = leg_max_x +2, 
                           xmax = leg_max_x +2.5, 
                           ymin = c(leg_max_y-2, leg_max_y-4), 
                           ymax = c(leg_max_y-1.5, leg_max_y-3.5 ), 
                           label_value = c("1", as.character(nrow(data.plot))), 
                           hjust = 0, vjust = 0, size = 2.5)
  text_data <- bind_rows(text_data, rank_title_data, arrow_text)
  
  # CREATE LEGEND for circle scores
  # circle legend
  if(F){
    cir_minimum_x <- x_min_score
    
    cir_legend_size <- 1
    cir_legend_space <- .1
    
    cir_legend_dat <-
      data.frame(
        value = seq(0, 1, by = .2),
        r = row_height/2*seq(0, 1, by = .2)
      )
    cir_legend_dat$r <- rescale(cir_legend_dat$r, to = c(0.05, 0.55), from = range(cir_legend_dat$r, na.rm = T))
    
    x0 <- vector("integer", nrow(cir_legend_dat))
    for(i in 1:length(x0)){
      if(i == 1){
        x0[i] <- cir_minimum_x + cir_legend_space + cir_legend_dat$r[i]
      }
      else {
        x0[i] <- x0[i-1] + cir_legend_dat$r[i-1] + cir_legend_space + cir_legend_dat$r[i]
      }
    }
    
    cir_legend_dat$x0 <- x0
    cir_legend_min_y <- leg_max_y-4
    cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r
    
    cir_legend_dat$colors <- NULL
    cir_maximum_x <- max(cir_legend_dat$x0)
    
    cir_title_data <- data_frame(xmin = cir_minimum_x,
                                 xmax = cir_maximum_x,
                                 ymin = leg_max_y -1,
                                 ymax = leg_max_y,
                                 label_value = "Score",
                                 hjust = 0, vjust = 0, fontface = "bold")
    
    cir_value_data <- data.frame(xmin = cir_legend_dat$x0 - cir_legend_dat$r,
                                 xmax = cir_legend_dat$x0 + cir_legend_dat$r,
                                 ymin = cir_legend_min_y,
                                 ymax = cir_legend_min_y +3,
                                 hjust = .5, vjust = 0, size = 2.5,
                                 label_value = ifelse(cir_legend_dat$value %in% c(0, 1),
                                                      paste0(cir_legend_dat$value*100, "%"), ""))
    
    circle_data <- bind_rows(circle_data, cir_legend_dat)
    text_data <- bind_rows(text_data, cir_title_data, cir_value_data)
    
    
    
  }
  
  minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE))
  
  # add to list
  plot_data[['row_info']] <- row_info
  plot_data[['column_info']] <- column_info
  plot_data[['palettes']] <- palettes
  plot_data[['row_pos']] <- row_pos
  plot_data[['column_pos']] <- column_pos
  plot_data[['rect_data']] <- rect_data
  plot_data[['text_data']] <- text_data
  plot_data[['image_data']] <- image_data
  plot_data[['arrow_data']] <- arrow_data
  
  print("Return a list object containing plot_data.")
  
  return(plot_data)
}


# function for plotting
scib_NicePlot <- function(plot_data, label_top3 = T){
  
  # initialization
  row_height <- 1.1
  row_space <- .1
  row_bigspace <- .5
  col_width <- 1.1
  col_space <- .2
  col_bigspace <- .5
  
  prepared_data_list <- plot_data
  
  g <-
    ggplot() +
    coord_equal(expand = FALSE) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_size_identity() +
    scale_linetype_identity() +
    cowplot::theme_nothing()
  
  # PLOT ROW BACKGROUNDS
  row_pos <- prepared_data_list[["row_pos"]]
  column_pos <- prepared_data_list[["column_pos"]]
  
  df <- row_pos %>% filter(colour_background)
  if (nrow(df) > 0) {
    g <- g + geom_rect(aes(xmin = min(column_pos$xmin)-.25, xmax = max(column_pos$xmax)+.25, ymin = ymin - (row_space / 2), ymax = ymax + (row_space / 2)), df, fill = "#DDDDDD")
  } 
  
  # PLOT CIRCLES
  if ( !is.null(prepared_data_list[["circle_data"]]) ) {
    g <- g + ggforce::geom_circle(aes(x0 = x0, y0 = y0, fill= colors, r = r), circle_data, size=.25)
  }
  
  # PLOT RECTANGLES
  if ( !is.null(prepared_data_list[["rect_data"]]) ) {
    rect_data <- prepared_data_list[["rect_data"]]
    # add defaults for optional values
    rect_data <- rect_data %>%
      add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
      mutate(border_colour = ifelse(border, border_colour, NA))
    
    g <- g + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colors, colour = border_colour, alpha = alpha), rect_data, size = .25)
  }
  
  # PLOT TEXT
  if ( !is.null(prepared_data_list[["text_data"]]) ) {
    text_data <- prepared_data_list[["text_data"]]
    # add defaults for optional values
    text_data <- text_data %>% 
      add_column_if_missing(
        hjust = .5,
        vjust = .5,
        size = 3,
        fontface = "plain",
        colors = "black",
        lineheight = 1,
        angle = 0
      ) %>%
      mutate(
        angle2 = angle / 360 * 2 * pi,
        cosa = cos(angle2) %>% round(2),
        sina = sin(angle2) %>% round(2),
        alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) + ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
        alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) + ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
        x = (1 - alphax) * xmin + alphax * xmax,
        y = (1 - alphay) * ymin + alphay * ymax
      ) %>%
      filter(label_value != "")
    # Set fontface for legend bold
    text_data[text_data$label_value == "Ranking", "fontface"] <- "bold"
    
    # subset text_data to left-aligned rows
    text_data_left <- text_data[which(text_data$group == "Methods" | text_data$group == "top3"), ]
    text_data <- text_data[-which(text_data$group == "Methods" | text_data$group == "top3"), ]
    
    g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = hjust, vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data)
    
    text_data_left[text_data_left$group == "Methods", "x"] <- text_data_left[text_data_left$group == "Methods", "x"] - 3
    if(label_top3){
      text_data_left[text_data_left$group == "top3", "x"] <- text_data_left[text_data_left$group == "top3", "xmin"] + .3
      text_data_left[text_data_left$group == "Methods", "x"] <- text_data_left[text_data_left$group == "Methods", "x"] + .5
    }
    g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = "left", vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data_left)
  }
  
  # PLOT ARROW RANKING
  if ( !is.null(prepared_data_list[["arrow_data"]]) ) {
    arrow_data <- prepared_data_list[["arrow_data"]]
    # add defaults for optional values
    arrow_data <- arrow_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
    
    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), arrow_data, arrow = arrow(length = unit(0.1, "cm")), lineend = "round", linejoin = "bevel")
  }
  
  # PLOT IMAGES
  # image_data$x = 1.2
  if( !is.null(prepared_data_list[["image_data"]]) ){
    image_data <- prepared_data_list[["image_data"]]
    for(r in 1:nrow(image_data)){
      img <- image_data[r, "image"]
      img_x <-  image_data[r, "x"]
      img_y <- image_data[r, "y"]
      g <- g+ cowplot::draw_image(image = img, x = img_x-.5, y = img_y-.5)
    }
    
  }
  
  # ADD SIZE
  # reserve a bit more room for text that wants to go outside the frame
  suppressWarnings({
    minimum_x <- min(column_pos$xmin, text_data$xmin, na.rm = TRUE)
    maximum_x <- max(column_pos$xmax, text_data$xmax, na.rm = TRUE)
    minimum_y <- min(row_pos$ymin, text_data$ymin, na.rm = TRUE)
    maximum_y <- max(row_pos$ymax, text_data$ymax, na.rm = TRUE)
  })
  
  minimum_x <- minimum_x - 0.5
  maximum_x <- maximum_x + 2
  minimum_y <- minimum_y - 1
  maximum_y <- maximum_y + 1.5
  
  g$width <- maximum_x - minimum_x
  g$height <- maximum_y - minimum_y
  
  g <- g + expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y))
  
  return(g)
}

# function for summary_tab to plot (OneClick)
scib_OneClick <- function(scib_summary_tab, label_top3 = T, Image_path = "./images"){
  preprocessed_data <- scib_score2tab(scib_summary_tab)
  prepared_data_list <- scib_tab2plot(preprocessed_data, label_top3 = T, Image_path = "./images")
  plot <- scib_NicePlot(plot_data = prepared_data_list, label_top3 = T)
  return(plot)
}


# functions to be used in the upper functions
add_column_if_missing <- function(df, ...) {
  column_values <- list(...)
  for (column_name in names(column_values)) {
    default_val <- rep(column_values[[column_name]], nrow(df))
    
    if (column_name %in% colnames(df)) {
      df[[column_name]] <- ifelse(is.na(df[[column_name]]), default_val, df[[column_name]])
    } else {
      df[[column_name]] <- default_val
    }
  }
  df
}


check_column_value <- function(data = data, 
                               column_name = "Output", 
                               allowed_values = c("Embedding", "Features", "Graph")
) {
  # 检查数据框中指定列的值是否都在允许的范围内
  if (!all(data[[column_name]] %in% allowed_values)) {
    # 找出不在允许范围内的值
    invalid_values <- data[[column_name]][!data[[column_name]] %in% allowed_values]
    # 报错并提示不在允许范围内的值
    stop("The following values in the '", column_name, "' column are not allowed: ",
         paste(invalid_values, collapse = ", "))
  }else{
    print(paste("Values in", column_name,"are corrected.", sep = " "))
  }
}

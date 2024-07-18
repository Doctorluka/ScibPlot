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

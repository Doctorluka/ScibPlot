# function for summary_tab to plot (OneClick)
scib_OneClick <- function(scib_summary_tab, label_top3 = T, Image_path = "./images"){
  preprocessed_data <- scib_score2tab(scib_summary_tab)
  prepared_data_list <- scib_tab2plot(preprocessed_data, label_top3 = T, Image_path = "./images")
  plot <- scib_NicePlot(plot_data = prepared_data_list, label_top3 = T)
  return(plot)
}


# data preprocessing from tab
scib_score2tab <- function(
    scib_summary_tab,
    Methods = "Methods",
    Features = "Features",
    Scaling = "Scaling",
    Output = "Output",
    Overall_Score = "Overall_Score",
    Bio_Score = "Bio_conservation",
    Batch_Score = "Batch_Correction",
    Image_path = "./images"
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

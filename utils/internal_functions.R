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

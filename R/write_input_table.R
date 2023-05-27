#' Write Input Table Helper Function
#'
#' This function generates table of input data
#' @param folder_name Character vector indicating folder or directory name to be used when outputting table
#' @param data Vector of observed values
#' @export
write_input_table <- function(folder_name, data){
  
  input_df <- 1:length(data) %>% as.data.frame()
  colnames(input_df) <- 'x'
  input_df$data <- data
  
  write.table(input_df, paste0(folder_name,'/input_data.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}
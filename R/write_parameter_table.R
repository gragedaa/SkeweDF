#' Write Parameter Table Helper Function
#'
#' This function generates table of optimized parameters
#' @param parameter_df Data frame of optimized parameters and other model function values (p0, Psi, etc)
#' @param folder_name Character vector indicating folder or directory name to be used when outputting table
#' @param model_fn_name Character vector used to indicate name of model function used for optimization
#' @param RGHD_m Int indicating m value of 2m-RGHD function if applicable
#' @export
write_parameter_table <- function(parameter_df, folder_name, model_fn_name, RGHD_m = 0){
  
  if(!(model_fn_name == 'RGHD')){
    fn_name <- str_replace(model_fn_name,'_',' ')
    
  }else{
    fn_name <- paste0('2m-RGHD (m=', RGHD_m, ')')
  }
  
  dir.create(paste0(folder_name,'/',fn_name))
  write.table(parameter_df, paste0(folder_name,'/',fn_name,'/',fn_name,' optimized_parameters.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}

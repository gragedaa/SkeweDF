#' Write Summary Table Helper Function
#'
#' This function generates summary statistics table of optimized parameters
#' @param parameter_df Data frame of optimized parameters and other model function values (p0, Psi, etc)
#' @param folder_name Character vector indicating folder or directory name to be used when outputting table
#' @param model_fn_name Character vector used to indicate name of model function used for optimization
#' @param RGHD_m Int indicating m value of 2m-RGHD function if applicable
#' @export
write_summary_table <- function(parameter_df, folder_name, model_fn_name, RGHD_m = 0){
  
  if(!(model_fn_name == 'RGHD')){
    fn_name <- str_replace(model_fn_name,'_',' ')
    
  }else{
    fn_name <- paste0('2m-RGHD (m=', RGHD_m, ')')
  }
  
  dir.create(paste0(folder_name,'/',fn_name))
  
  summary_df <- colnames(parameter_df) %>% as.data.frame()
  colnames(summary_df) <- 'Variable'
  summary_df$n <- nrow(parameter_df)
  summary_df$mean <- colMeans(parameter_df)
  summary_df$median <- colMedians(parameter_df %>% as.matrix())
  summary_df$std <- colSds(parameter_df %>% as.matrix())
  summary_df$mean_error_alpha_0.05 <- qt(1-0.05, summary_df$n-1) * summary_df$std / sqrt(summary_df$n)
  summary_df$mean_CI_alpha_0.05 <- paste0((summary_df$mean - summary_df$mean_error_alpha_0.05) %>% round(5), ', ',(summary_df$mean + summary_df$mean_error_alpha_0.05) %>% round(5))
  summary_df$median_CI_lower_rank <- ((summary_df$n/2) - (1.96 * sqrt(summary_df$n) / 2)) %>% round(0)
  summary_df$median_CI_upper_rank <- (1+ (summary_df$n/2) + (1.96 * sqrt(summary_df$n) / 2)) %>% round(0)
  summary_df$median_CI_approx_0.95 <- ''
  for(i in 1:length(parameter_df)){
    tmp <- sort(parameter_df[i] %>% unlist())
    summary_df$median_CI_approx_0.95[i] <- paste0(tmp[ summary_df$median_CI_lower_rank[i]]%>%round(5),', ',tmp[ summary_df$median_CI_upper_rank[i]]%>%round(5))
  }
  
  write.table(summary_df, paste0(folder_name,'/',fn_name,'/',fn_name,' parameter_summary.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}

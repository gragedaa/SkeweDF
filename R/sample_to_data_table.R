#' Sample to Data Table
#' 
#' Helper function which converts a sample vector to a frequency dataframe.
#' @param sample Vector
#'
#' @export
sample_to_data_table <- function(sample){
  out <- data.frame(table(sample))
  colnames(out) <- c('x','Freq')
  out[,1] <- out[,1] %>% as.character() %>% as.numeric()
  tmp <- 1:max(out[,1]) %>% as.data.frame()
  colnames(tmp) <- 'x'
  out <- left_join(tmp,out)
  out$Freq[is.na(out$Freq)] <- 0
  return(out)
}
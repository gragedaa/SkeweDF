#' Data Table to Sample
#' 
#' Helper function which converts a frequency dataframe to a sample vector
#' @param dt DataFrame. Frequency dataframe 
#'
#' @export
data_table_to_sample <- function(dt){
  out <- c()
  for(i in 1:nrow(dt)){
    out <- c(out,rep(dt[i,1], dt[i,2]))
  }
  return(out)
}
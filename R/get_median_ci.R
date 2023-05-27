#' Get Median Confidence Interval Function
#'
#' This function generates a vector of ranked 95% confidence interval based on median of data.
#' @param data Data to get confidence interval from
#' @export
get_median_CI <- function(data){
  n <- length(data)
  lower_rank <- (n/2) - (1.96 * sqrt(n) / 2)
  upper_rank <- 1 + (n/2) + (1.96 * sqrt(n) / 2)
  lower_rank <- lower_rank %>% round(0)
  upper_rank <- upper_rank %>% round(0)
  
  data <- sort(data)
  return(c(data[lower_rank], data[upper_rank]))
}
#' Get Mean Confidence Interval Function
#'
#' This function generates a vector of confidence interval based on mean of data.
#' @param data Data to get confidence interval from
#' @param alpha Alpha for confidence interval calculation
#' @export
get_CI <- function(data, alpha){
  std <- sd(data)
  mean <- mean(data)
  sample_size <- length(data)
  mean_error <- qt(1-alpha, sample_size-1) * std / sqrt(sample_size)
  
  return(c(mean - mean_error, mean + mean_error))
}
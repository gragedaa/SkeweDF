#' Right-Tail Cumulative Distribution Function
#'
#' This function generates a vector of the right-tail cumulative distribution function of a given vector of values.
#' @param x Length of vector to be generated.
#' @examples
#' x <- c(1,2,3,4,5)
#' right_tail_cdf(x)
#' @export
right_tail_cdf <- function(x){
  output_cdf <- x[length(x):1]
  #output_cdf <- output_cdf / sum(output_cdf)
  for(i in 2:length(output_cdf)){
    output_cdf[i] <- output_cdf[i] + output_cdf[i-1]
  }
  output_cdf <- output_cdf[length(output_cdf):1]
}
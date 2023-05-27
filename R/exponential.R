#' Exponential Distribution Function
#'
#' This function generates a vector of n length of the Exponential distribution with parameters a and b.
#' @param n Length of vector to be generated.
#' @param a Parameter of the Exponential distribution function
#' @param b Parameter of the Exponential distribution function
#' @examples
#' Exponential(100, 10000, 0.8)
#' @export
Exponential <- function(n, a , b){
  p <- a * exp(-b * (1:n))
  return(p)
}
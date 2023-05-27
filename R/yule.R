#' Yule Distribution Function
#'
#' This function generates a vector of n length of the Yule distribution with parameter rho.
#' @param n Length of vector to be generated.
#' @param rho Parameter of the Yule distribution function
#' @examples
#' Yule(100, 3)
#' @export
Yule <- function(n, rho){
  p <- beta(1:n,rho+1)
  p <- p * rho
  p <- p / sum(p)
  return(p)
}
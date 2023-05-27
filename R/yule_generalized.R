#' Generalized Yule Distribution Function
#'
#' This function generates a vector of n length of the Generalized Yule distribution with parameters rho and alpha.
#' @param n Length of vector to be generated.
#' @param rho Parameter of the Generalized Yule distribution function
#' @param alpha Parameter of the Generalized Yule distribution function: 0 <= alpha < 1
#' @examples
#' Generalized_Yule(100, 3, 0.1)
#' @export
Generalized_Yule <- function(n, rho, alpha){
  out <- rho / (1 - alpha^rho)
  two <- Ibeta(alpha,1:n, rho+1)
  out <- out * two
  return(out / sum(out))
}
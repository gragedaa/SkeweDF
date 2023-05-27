#' Lorentzian Distribution Function calculation
#'
#' This function calculates value of Lorentzian function at x
#' @param x Index of function
#' @param gamma Parameter of the Lorenzian distribution function
#' @param x0 Parameter of the Lorenzian distribution function indicating center of function
#' @param c Parameter of the Lorenzian distribution function indicating center of function
#' @examples
#' Lorentzian_calc(5, 5.5, 6, 2)
#' @export
Lorentzian_calc <- function(x, gamma, x0, c){
  
  a <- abs(x- x0) ^ c
  
  b <- gamma / 2
  
  output <- b / (a + b^c)
  return(output)
}

#' Lorentzian Distribution Function
#'
#' This function generates a vector of n length of the Lorentzian distribution
#' @param n Length of vector to be generated.
#' @param gamma Parameter of the Lorenzian distribution function
#' @param x0 Parameter of the Lorenzian distribution function indicating center of function
#' @param c Parameter of the Lorenzian distribution function indicating center of function
#' @examples
#' Lorentzian_calc(5, 5.5, 6, 2)
#' @export
Lorentzian <- function(n, gamma, x0, c){
  out <- lapply(1:n, function(i){
    return(Lorentzian_calc(i, gamma, x0, c))
  })
  out <- unlist(out)
  out <- out / sum(out)
  return(out)
}
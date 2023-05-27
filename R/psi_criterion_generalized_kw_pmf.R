#' Psi Criterion Generalized KW PMF
#'
#' This function calculates the Psi criterion goodness of cit metric given a set of parameters for the probability mass function of the Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @export
psi_criterion_generalized_kw_pmf <- function(x, d, data, left_trunc, right_trunc){
  a <- x[1:d]
  b <- x[(d+1):(length(x)-1)]
  theta <- x[length(x)]
  
  model <- Kolmogorov_Waring(right_trunc, a, b, theta)
  model <- model[-1]
  model <- model[left_trunc:right_trunc]
  model <- model / sum(model)
  model <- model * sum(data)
  
  
  return(psi_criterion(data, model, length(x)) * -1)
}
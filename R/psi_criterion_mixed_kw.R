#' Psi Criterion Mixed KW
#'
#' This function calculates the Psi criterion goodness of fit metric given a set of parameters for a mixed KW function
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @export
psi_criterion_mixed_kw <- function(x, d, data, left_trunc, right_trunc, pmf_weight, weighted_rt){
  pmf <- psi_criterion_mixed_kw_pmf(x, d, data, left_trunc, right_trunc)
  cdf <- psi_criterion_mixed_kw_cdf(x, d, data, left_trunc, right_trunc, weighted_rt)
  
  return((pmf_weight * pmf) + ((1-pmf_weight)*cdf))
}
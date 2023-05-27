#' Psi Criterion Generalized KW
#'
#' This function calculates the Psi criterion goodness of fit metric given a set of parameters for the Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param pmf_weight Numeric of weight given to probability mass function for generation of Psi Criterion. For example, if pmf_weight <- 0.5, 50 percent of the Psi Criterion value will be attributed to the probability mass function while the other 50 percent will be attributed to the right-tail cumulative distribution function.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param ratio Boolean used to determine if parameters provide are a/b ratios rather than absolute values.
#' @export
psi_criterion_generalized_kw <- function(x, d, data, left_trunc, right_trunc, pmf_weight, weighted_rt, ratio = FALSE){
  if(ratio){
    x[1:d] <- x[1:d] * x[(d+1):(d+d)]
    print(x)
  }
  
  pmf <- psi_criterion_generalized_kw_pmf(x, d, data, left_trunc, right_trunc)
  cdf <- psi_criterion_generalized_kw_cdf(x, d, data, left_trunc, right_trunc, weighted_rt)
  
  return((pmf_weight * pmf) + ((1-pmf_weight)*cdf))
}

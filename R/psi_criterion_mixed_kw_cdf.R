#' Psi Criterion Mixed KW CDF
#'
#' This function calculates the Psi criterion goodness of fit metric given a set of parameters for the cumulative distribution function of the a mixed Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @export
psi_criterion_mixed_kw_cdf <- function(x, d, data, left_trunc, right_trunc, weighted_rt){
  a1 <- x[1:d[1]]
  b1 <- x[(d[1]+1):(2*d[1])]
  theta1 <- x[(2*d[1]+1)]
  
  a2 <- x[(2*d[1]+2):(2*d[1]+2+d[2]-1)]
  b2 <- x[(2*d[1]+2+d[2]):(2*d[1]+2+2*d[2]-1)]
  theta2 <- x[2*d[1]+2+2*d[2]]
  
  model1_weight <- unlist(x[2*d[1]+2+2*d[2]+1])
  
  model1 <- Kolmogorov_Waring(right_trunc, unlist(a1), unlist(b1), unlist(theta1))
  model2 <- Kolmogorov_Waring(right_trunc, unlist(a2), unlist(b2), unlist(theta2))
  model <- additive_mixed_model(model1, model2, model1_weight)
  model <- model[-1]
  model <- model[left_trunc:right_trunc]
  model <- model / sum(model)
  model <- model * sum(data)
  
  if(weighted_rt){
    data <- weighted_right_tail_cdf(data)
    model <- weighted_right_tail_cdf(model)
  }
  else{
    data <- right_tail_cdf(data)
    model <- right_tail_cdf(model)
  }
  
  
  return(psi_criterion(data, model, length(x)) * -1)
}

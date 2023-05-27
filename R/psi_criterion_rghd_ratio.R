#' Psi Criterion for RGHD parameter ratios
#'
#' This function generates the Psi Criterion goodness of fit value given an empirical distribution for the 2m-RGHD function. Parameters r and q/r ratios are given, as well as desired weight of pmf and use of the weighted right-tail cumulative distribution function.
#' @param params Vector of parameters for model_fn, not including n. For example, for 2m-RGHD (m=2), params <- c(3, 5, 0.3, 1.5). In this case r1 = 3, r2 = 5, q1/r1 = 0.3, and q2/r2 = 1.5
#' @param data Vector of observed values
#' @param m m parameter for 2m-RGHD function
#' @param pmf_weight Numeric of weight given to probability mass function for generation of Psi Criterion. For example, if pmf_weight <- 0.5, 50 percent of the Psi Criterion value will be attributed to the probability mass function while the other 50 percent will be attributed to the right-tail cumulative distribution function.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @examples
#' obs_data <- c(100,75,20,1)
#' parameters <- c(3, 5, 0.3, 1.5)
#' psi <- psi_criterion_RGHD_ratio(parameters, obs_data, 2)
#' @export
psi_criterion_RGHD_ratio <- function(params, data, m, pmf_weight = 0, weighted_rt = FALSE, left_trunc = 1, right_trunc = left_trunc + length(data) - 1){ ## multistart can't maximize atm so I just negate the criterion for simple minimization
  
  model <- RGHD(right_trunc, m, unlist(params[1:m]), unlist(params[1:m]) * unlist(params[(m+1):(m+m)]))
  model <- model[left_trunc:right_trunc]
  model <- model / sum(model)
  
  if(weighted_rt){
    data_cdf <- weighted_right_tail_cdf(data)
    model_cdf <- weighted_right_tail_cdf(model)
  }
  else{
    data_cdf <- right_tail_cdf(data)
    model_cdf <- right_tail_cdf(model)
  }
  
  return((((pmf_weight)*psi_criterion(data, model, m*2)) + ((1-pmf_weight)*psi_criterion(data_cdf, model_cdf, m*2))) * -1)
}
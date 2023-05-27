#' Psi Criterion given a function
#'
#' This function generates the Psi Criterion goodness of fit value given an empirical distribution. The function and parameters are given, as well as desired weight of pmf and use of the weighted right-tail cumulative distribution function.
#' @param params Vector of parameters for model_fn, not including n. For example, for Generalized_Yule(n, rho, alpha), params will be c(rho, alpha)
#' @param data Vector of observed values
#' @param model_fn Function of theoretical model to be used. For example, for Generalized_Yule(n, rho, alpha), model_fn <- Generalied_Yule
#' @param pmf_weight Numeric of weight given to probability mass function for generation of Psi Criterion. For example, if pmf_weight <- 0.5, 50 percent of the Psi Criterion value will be attributed to the probability mass function while the other 50 percent will be attributed to the right-tail cumulative distribution function.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @examples
#' obs_data <- c(100,75,20,1)
#' parameters <- c(1,2,0.8)
#' psi <- psi_criterion_function(parameters, obs_data, Kolmogorov_Waring)
#' @export
psi_criterion_function <- function(params, data, model_fn, pmf_weight = 0, weighted_rt = FALSE, left_trunc = 1, right_trunc = left_trunc + length(data) - 1){
  
  model_function <- match.fun(model_fn)
  
  model <- invoke(model_function, c(right_trunc, params) %>% unlist() %>% unname())
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
  
  return((((pmf_weight)*psi_criterion(data, model, length(params))) + ((1-pmf_weight)*psi_criterion(data_cdf, model_cdf, length(params)))) * -1)
}
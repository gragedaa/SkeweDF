#' Psi Criterion
#'
#' This function generates the Psi Criterion goodness of fit value given an empirical distribution, theoretical modeled distribution, and number of parameters in the theoretical distribution.
#' @param data Vector of observed values
#' @param model Vector of theoretical values to be compared
#' @param n_parameters Number of parameters of function used to generate model
#' @examples
#' obs_data <- c(100,75,20,1)
#' model_data <- Kolmogorov_Waring(length(obs_data), 2, 3, 0.9)
#' psi <- psi_criterion(obs_data, model_data, 3)
#' @export
psi_criterion <- function(data, model, n_parameters){
  var_data <- data - mean(data)
  var_data <- var_data * var_data
  
  diff <- data - model
  diff <- diff * diff
  
  psi <- log(sum(var_data) / sum(diff)) - (2 * n_parameters / length(data))
  
  return(psi)
}
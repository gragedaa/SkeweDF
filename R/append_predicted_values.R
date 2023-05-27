#' Append Predicted Values
#'
#' This function will take in a left truncated data vector and append missing values according to the given model
#' @param data Vector of observed values
#' @param model Vector of optimized model
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @export
append_predicted_values <- function(data, model, left_trunc = 1, right_trunc = left_trunc + length(data) - 1){
  predict <- model[1:left_trunc] / sum(model[-(1:left_trunc)]) * sum(data)
  full_data <- c(predict, data)
  return(full_data)
}
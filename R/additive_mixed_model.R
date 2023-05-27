#' Additive Mixed Model
#'
#' This function will return the weighted sum between two models of the same length
#' @param model1 Vector.
#' @param model2 Vector.
#' @param model1_weight numeric. Weight of model1
#' 
#' @export
additive_mixed_model <- function(model1, model2, model1_weight){
  output <- (model1_weight * model1) + ((1-model1_weight) * model2)
  return(output)
}
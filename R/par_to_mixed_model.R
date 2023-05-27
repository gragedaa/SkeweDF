#' Parameters to Mixed Model
#'
#' Helper function which inputs a vector of parameters output from fit_mixed_kw and creates a model
#' @param n Int. Length of model vector
#' @param x Vector. Parameters
#' @param d Vector. Dimensions of mixed model
#'
#' @export
par_to_mixed_model <- function(n, x, d){
  a1 <- x[1:d[1]]
  b1 <- x[(d[1]+1):(2*d[1])]
  theta1 <- x[(2*d[1]+1)]
  
  a2 <- x[(2*d[1]+2):(2*d[1]+2+d[2]-1)]
  b2 <- x[(2*d[1]+2+d[2]):(2*d[1]+2+2*d[2]-1)]
  theta2 <- x[2*d[1]+2+2*d[2]]
  
  model1_weight <- unlist(x[2*d[1]+2+2*d[2]+1])
  
  model1 <- Kolmogorov_Waring(n, unlist(a1), unlist(b1), unlist(theta1))
  model2 <- Kolmogorov_Waring(n, unlist(a2), unlist(b2), unlist(theta2))
  model <- additive_mixed_model(model1, model2, model1_weight)
  
  return(model)
}

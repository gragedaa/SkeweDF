#' Psi Criterion for RGHD parameter ratios
#'
#' This function generates the Psi Criterion goodness of fit value given an empirical distribution for the 2m-RGHD function. Parameters r and q/r ratios are given, as well as desired weight of pmf and use of the weighted right-tail cumulative distribution function.
#' @param params Vector of parameter for the model function
#' @param model_fn_name name of function as a character vector
#' @examples
#' params <- c(2, 3, 0.9)
#' get_p0(params, 'Kolmogorov Waring')
#' @export
get_p0 <- function(params, model_fn_name){
  if(model_fn_name == 'Yule'){
    return(NA)
  }else if(model_fn_name == 'Generalized_Yule'){
    return(NA)
  }else if(model_fn_name == 'Kolmogorov_Waring'){
    return(Kolmogorov_Waring_P0(params[1],params[2],params[3]))
  }else if(model_fn_name == 'RGHD'){
    m <- length(params)/2
    return(RGHD_P0_calc(100,m,params[1:m] %>% unlist(), params[(m+1):(m+m)] %>% unlist()))
  }
}
#' Fit Generalized KW
#'
#' This function will take in a set of parameters and calculate goodness of fit based on the Psi criterion metric
#' @param param_bounds A list of sequences which indicate space where parameters should be generated and fit
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param pmf_weight Float used to determine weight of probability mass function for goodness of fit calculation relative to the cumulative distribution function
#' @param par_chunk Integer used to indicate number of optimization chunks to be run. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param par_chunk_size Integer used to indicate number of starting parameters to be generated and optimized in a given chunk. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param n_cores Integer used to indicate number of cores to be used for this function if a socket cluster object is not defined.
#' @param clust socket cluster object from 'parallel::makeCluster()'. This is used if you have already generated a socket cluster object and would like to run this functoin on it. If no object is defined, one will be made for this function call.
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param ratio Boolean used to determine if parameters provide are a/b ratios rather than absolute values.
#' @export
fit_generalized_kw <- function(param_bounds, d, data, weighted_rt = FALSE, pmf_weight = 0.0, par_chunk = 100, par_chunk_size = 10, n_cores = 1, clust, left_trunc = 1, right_trunc = left_trunc+length(data)-1, ratio = FALSE){
  
  data <- data / sum(data)
  
  standalone <- FALSE
  
  
  if(missing(clust)){
    standalone <- TRUE
    clust <- makeCluster(n_cores)
    clusterExport(clust, varlist = c('data','param_bounds','weighted_rt', 'par_chunk_size','d','left_trunc','right_trunc', 'pmf_weight','ratio'), envir = environment())
    clusterCall(clust, function() library(SkeweDF))
  }
  
  clusterExport(clust, varlist = c('data','param_bounds','weighted_rt', 'par_chunk_size','d','left_trunc','right_trunc', 'pmf_weight','ratio'), envir = environment())
  
  parameters <- parLapply(clust,1:par_chunk, function(q){
    par_mat <- lapply(1:par_chunk_size, function(i){
      par_vec <- c(rep(0, length(param_bounds)))
      for(i in 1:length(param_bounds)){
        par_vec[i] <- sample(param_bounds[[i]], size = 1, replace = TRUE)
      }
      return(par_vec)
    }) %>% as.data.frame() %>% t() %>% as.matrix()
    
    param_lower <- c(rep(0, length(param_bounds)))
    param_upper <- c(rep(0, length(param_bounds)))
    
    for(i in 1:length(param_bounds)){
      param_lower[i] <- min(param_bounds[[i]])
      param_upper[i] <- max(param_bounds[[i]])
    }
    
    fn_parameters <- multistart(parmat = par_mat,fn = psi_criterion_generalized_kw,method = 'L-BFGS-B',
                                lower = param_lower, upper = param_upper,
                                data = data, d = d, pmf_weight = pmf_weight, weighted_rt = weighted_rt,
                                left_trunc = left_trunc, right_trunc = right_trunc, ratio = ratio);
    
    fn_parameters$fevals <- NULL
    fn_parameters$gevals <- NULL
    fn_parameters$convergence <- NULL
    
    for(i in 1:d){
      if(ratio){
        colnames(fn_parameters)[i] <- paste0('ab_ratio',i)
      }else{
        colnames(fn_parameters)[i] <- paste0('a',i)
      }
    }
    for(i in (d+1):length(fn_parameters)){
      colnames(fn_parameters)[i] <- paste0('b',i-d)
    }
    
    #fn_parameters[1:m] <- t(apply(fn_parameters[1:m], 1, FUN=function(x) sort(x)))
    #fn_parameters[(m+1):(length(fn_parameters)-2)] <- t(apply(fn_parameters[(m+1):(length(fn_parameters)-2)], 1, FUN=function(x) sort(x)))
    
    colnames(fn_parameters)[length(fn_parameters)-1] <- 'theta'
    colnames(fn_parameters)[length(fn_parameters)] <- 'Psi'
    fn_parameters$Psi <- fn_parameters$Psi * -1
    fn_parameters$pmf_weight <- pmf_weight
    return(fn_parameters)
  }) %>% bind_rows()
  
  print('Complete')
  
  if(standalone){
    stopCluster(clust)
  }
  
  parameters <- arrange(parameters, desc(parameters$Psi))
  if(!ratio){
    if(length(1:d) > 1){
      parameters[1:d] <- t(apply(parameters[1:d], 1, FUN=function(x) sort(x)))
    }
    if(length((d+1):(length(parameters)-3)) > 1){
      parameters[(d+1):(length(parameters)-3)] <- t(apply(parameters[(d+1):(length(parameters)-3)], 1, FUN=function(x) sort(x)))
    }
  }
  
  
  #parameters <- parameter_post_processing(parameters, 'Generalized KW', data)
  
  return(parameters)
  
}

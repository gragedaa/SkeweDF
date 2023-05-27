#' Local optimization of the 2m-RGHD function given empirical data, r bounds, and q/r bounds.
#'
#' This function generates a table of optimized parameters and Psi Criterion for a given function within specified starting parameter bounds. This function uses Limited Memory BFGS as it's gradient descent algorithm.
#' @param param_bounds A list of sequences which indicate space where parameters should be generated and fit
#' @param data Vector of observed values
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param par_chunk Integer used to indicate number of optimization chunks to be run. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param par_chunk_size Integer used to indicate number of starting parameters to be generated and optimized in a given chunk. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param n_cores Integer used to indicate number of cores to be used for this function if a socket cluster object is not defined.
#' @param clust socket cluster object from 'parallel::makeCluster()'. This is used if you have already generated a socket cluster object and would like to run this functoin on it. If no object is defined, one will be made for this function call.
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @export
local_fit_RGHD_ratio <- function(param_bounds, data, weighted_rt = FALSE, pmf_weight = 0.0, par_chunk = 100, par_chunk_size = 10, n_cores = 1, clust, left_trunc = 1, right_trunc = left_trunc+length(data)-1){
  
  if(length(param_bounds) %% 2 != 0){
    stop('Must have even number of parameters!')
  }
  
  skeweDF_gVar_data <- data / sum(data)
  
  standalone <- FALSE
  
  skeweDF_gVar_m <- length(param_bounds) / 2
  skeweDF_gVar_param_bounds <- param_bounds
  skeweDF_gVar_weighted_rt <- weighted_rt
  skeweDF_gVar_left_trunc <- left_trunc
  skeweDF_gVar_right_trunc <- right_trunc
  skeweDF_gVar_par_chunk_size <- par_chunk_size
  
  if(missing(clust)){
    standalone <- TRUE
    clust <- makeCluster(n_cores)
    clusterExport(clust, varlist = c('skeweDF_gVar_data','skeweDF_gVar_param_bounds','skeweDF_gVar_weighted_rt', 'skeweDF_gVar_par_chunk_size','skeweDF_gVar_m','skeweDF_gVar_left_trunc','skeweDF_gVar_right_trunc', 'pmf_weight'), envir = environment())
    clusterCall(clust, function() library(SkeweDF))
  }
  
  
  clusterExport(clust, varlist = c('skeweDF_gVar_data','skeweDF_gVar_param_bounds','skeweDF_gVar_weighted_rt', 'skeweDF_gVar_par_chunk_size','skeweDF_gVar_m','skeweDF_gVar_left_trunc','skeweDF_gVar_right_trunc', 'pmf_weight'), envir = environment())
  
  parameters <- parLapply(clust,1:par_chunk, function(q){
    par_mat <- lapply(1:skeweDF_gVar_par_chunk_size, function(i){
      par_vec <- c(rep(0, length(skeweDF_gVar_param_bounds)))
      for(i in 1:length(skeweDF_gVar_param_bounds)){
        par_vec[i] <- sample(skeweDF_gVar_param_bounds[[i]], size = 1, replace = TRUE)
      }
      return(par_vec)
    }) %>% as.data.frame() %>% t() %>% as.matrix()
    
    param_lower <- c(rep(0, length(skeweDF_gVar_param_bounds)))
    param_upper <- c(rep(0, length(skeweDF_gVar_param_bounds)))
    
    for(i in 1:length(skeweDF_gVar_param_bounds)){
      param_lower[i] <- min(skeweDF_gVar_param_bounds[[i]])
      param_upper[i] <- max(skeweDF_gVar_param_bounds[[i]])
    }
    
    fn_parameters <- multistart(parmat = par_mat,fn = psi_criterion_RGHD_ratio,method = 'L-BFGS-B',
                                lower = param_lower, upper = param_upper,
                                data = skeweDF_gVar_data, m = skeweDF_gVar_m, pmf_weight = pmf_weight, weighted_rt = skeweDF_gVar_weighted_rt,
                                left_trunc = skeweDF_gVar_left_trunc, right_trunc = skeweDF_gVar_right_trunc);
    
    for(i in 1:skeweDF_gVar_m){
      colnames(fn_parameters)[i] <- paste0('r',i)
      colnames(fn_parameters)[i+skeweDF_gVar_m] <- paste0('q',i)
    }
    colnames(fn_parameters)[skeweDF_gVar_m+skeweDF_gVar_m+1] <- 'Psi_RTCDF'
    fn_parameters$Psi_RTCDF <- fn_parameters$Psi_RTCDF * -1
    return(fn_parameters)
  }) %>% bind_rows()
  
  print('Complete')
  
  if(standalone){
    stopCluster(clust)
  }
  
  parameters <- arrange(parameters, desc(parameters$Psi_RTCDF))
  
  for(i in 1:skeweDF_gVar_m){
    parameters[skeweDF_gVar_m+i] <- parameters[skeweDF_gVar_m+i] * parameters[i]
  }
  
  parameters$fevals <- NULL
  parameters$gevals <- NULL
  parameters$convergence <- NULL
  
  parameters <- parameter_post_processing(parameters, 'RGHD', skeweDF_gVar_data)
  
  return(parameters)
  
}

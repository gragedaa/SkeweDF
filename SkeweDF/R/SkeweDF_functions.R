## usethis namespace: start
#' @useDynLib SkeweDF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' Yule Distribution Function
#'
#' This function generates a vector of n length of the Yule distribution with parameter rho.
#' @param n Length of vector to be generated.
#' @param rho Parameter of the Yule distribution function
#' @examples
#' Yule(100, 3)
#' @export
Yule <- function(n, rho){
  p <- beta(1:n,rho+1)
  p <- p * rho 
  p <- p / sum(p)
  return(p)
}

#' Generalized Yule Distribution Function
#'
#' This function generates a vector of n length of the Generalized Yule distribution with parameters rho and alpha.
#' @param n Length of vector to be generated.
#' @param rho Parameter of the Generalized Yule distribution function
#' @param alpha Parameter of the Generalized Yule distribution function: 0 <= alpha < 1
#' @examples
#' Generalized_Yule(100, 3, 0.1)
#' @export
Generalized_Yule <- function(n, rho, alpha){
  out <- rho / (1 - alpha^rho)
  two <- Ibeta(alpha,1:n, rho+1)
  out <- out * two
  return(out / sum(out))
}

#' Weighted Right-Tail Cumulative Distribution Function
#'
#' This function generates a vector of the weighted right-tail cumulative distribution function of a given vector of values. The weight of of each variable is determined by its position in the vector. For example, with a vector of length 5, element 5 will have weight 5/(5+4+3+2+1). Element 1 will have weight 1/(5+4+3+2+1)
#' @param x Length of vector to be generated.
#' @examples
#' x <- c(1,2,3,4,5)
#' weighted_right_tail_cdf(x)
#' @export
weighted_right_tail_cdf <- function(x){
  output_cdf <- x[length(x):1]
  output_cdf <- output_cdf * length(x):1
  output_cdf <- output_cdf / sum(output_cdf)
  for(i in 2:length(output_cdf)){
    output_cdf[i] <- output_cdf[i] + output_cdf[i-1]
  }
  output_cdf <- output_cdf[length(output_cdf):1]
}

#' Right-Tail Cumulative Distribution Function
#'
#' This function generates a vector of the right-tail cumulative distribution function of a given vector of values.
#' @param x Length of vector to be generated.
#' @examples
#' x <- c(1,2,3,4,5)
#' right_tail_cdf(x)
#' @export
right_tail_cdf <- function(x){
  output_cdf <- x[length(x):1]
  #output_cdf <- output_cdf / sum(output_cdf)
  for(i in 2:length(output_cdf)){
    output_cdf[i] <- output_cdf[i] + output_cdf[i-1]
  }
  output_cdf <- output_cdf[length(output_cdf):1]
}

#' Psi Criterion
#'
#' This function generates the Psi Criterion goodness of fit value given an empirical distribution, theoretical modeled distribution, and number of parameters in the theoretical distribution.
#' @param data Vector of observed values 
#' @param model Vector of theoretical values to be compared
#' @param n_parameters Number of parameters of function used to generate model
#' @examples
#' obs_data <- c(100,75,20,1)
#' model_data <- theoretical_model(param1,param2,param3)
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


#' Psi Criterion given a function 
#'
#' This function generates the Psi Criterion goodness of fit value given an empirical distribution. The function and parameters are given, as well as desired weight of pmf and use of the weighted right-tail cumulative distribution function.
#' @param params Vector of parameters for model_fn, not including n. For example, for Generalized_Yule(n, rho, alpha), params will be c(rho, alpha)
#' @param data Vector of observed values 
#' @param model_fn Function of theoretical model to be used. For example, for Generalized_Yule(n, rho, alpha), model_fn <- Generalied_Yule
#' @param pmf_weight Numeric of weight given to probability mass function for generation of Psi Criterion. For example, if pmf_weight <- 0.5, 50% of the Psi Criterion value will be attributed to the probability mass function while the other 50% will be attributed to the right-tail cumulative distribution function.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @examples
#' obs_data <- c(100,75,20,1)
#' parameters <- c(1,2,0.8)
#' psi_criterion_function <- psi_criterion_function(parameters, obs_data, Kolmogorov_Waring, 0, F)
#' @export
psi_criterion_function <- function(params, data, model_fn, pmf_weight = 0, weighted_rt = F){
  
  model_function <- match.fun(model_fn)
  
  model <- invoke(model_function, c(length(data), params) %>% unlist() %>% unname())
  
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

#' Psi Criterion for RGHD parameter ratios
#'
#' This function generates the Psi Criterion goodness of fit value given an empirical distribution for the 2m-RGHD function. Parameters r and q/r ratios are given, as well as desired weight of pmf and use of the weighted right-tail cumulative distribution function.
#' @param params Vector of parameters for model_fn, not including n. For example, for 2m-RGHD (m=2), params <- c(3, 5, 0.3, 1.5). In this case r1 = 3, r2 = 5, q1/r1 = 0.3, and q2/r2 = 1.5
#' @param data Vector of observed values 
#' @param m m parameter for 2m-RGHD function
#' @param pmf_weight Numeric of weight given to probability mass function for generation of Psi Criterion. For example, if pmf_weight <- 0.5, 50% of the Psi Criterion value will be attributed to the probability mass function while the other 50% will be attributed to the right-tail cumulative distribution function.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @examples
#' obs_data <- c(100,75,20,1)
#' parameters <- c(3, 5, 0.3, 1.5)
#' psi_criterion_RGHD_ratio <- psi_criterion_function(parameters, obs_data, 2, 0, F)
#' @export
psi_criterion_RGHD_ratio <- function(params, data, m, pmf_weight, weighted_rt){ ## multistart can't maximize atm so I just negate the criterion for simple minimization
  
  model <- RGHD(length(data), m, unlist(params[1:m]), unlist(params[1:m]) * unlist(params[(m+1):(m+m)]))
  
  if(weighted_rt){
    data_cdf <- weighted_right_tail_cdf(data)
    model_cdf <- weighted_right_tail_cdf(model)
  }
  else{
    data_cdf <- right_tail_cdf(data)
    model_cdf <- right_tail_cdf(model)
  }
  
  return((((pmf_weight)*psi_criterion(data, model, m*2)) + ((1-pmf_weight)*psi_criterion(data_cdf, model_cdf, m*2))) * -1)
  #return(psi_criterion(data, model, m*2) * psi_criterion(data_cdf, model_cdf, m*2) * -1)
}

#' Local optimization of a given function given empirical data and parameter bounds
#'
#' This function generates a table of optimized parameters and Psi Criterion for a given function within specified starting parameter bounds. This function uses Limited Memory BFGS as it's gradient descent algorithm.
#' @param param_bounds A list of sequences which indicate space where parameters should be generated and fit
#' @param data Vector of observed values 
#' @param model_fn Function of theoretical model to be used. For example, for Generalized_Yule(n, rho, alpha), model_fn <- Generalied_Yule
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param par_chunk Integer used to indicate number of optimization chunks to be run. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param par_chunk_size Integer used to indicate number of starting parameters to be generated and optimized in a given chunk. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param n_cores Integer used to indicate number of cores to be used for this function if a socket cluster object is not defined.
#' @param clust socket cluster object from 'parallel::makeCluster()'. This is used if you have already generated a socket cluster object and would like to run this functoin on it. If no object is defined, one will be made for this function call.
#' @examples
#' obs_data <- c(16362,7910,4563,2880,1992,1242,936,734,553,346,253,296,185,210,120,121,119,94,87,55,47,60,67,26,24,54,32,31,16,19,19,20,11,22,21,25,8,8,28,9,4,4,2,8,6,5,1,9,0,6,2,4,2,1,0,4,3,2,6,5,3,1,1,0,3,4,0,1,0,0,0,1,1,0,0,0,0,0,0,3,0,0,0,0,1,0,0,2,1,5,0,0,0,1,1,1,0,0,0,0,0,0,0,5,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4)
#' param_bounds <- list(seq(0.1,5,0.1),seq(0.1,5,0.1),seq(0.01,1,0.01))
#' parameters <- local_fit_function(param_bounds = param_bounds, data = obs_data, model_fn = Kolmogorov_Waring, n_cores = 12)
#' 
#' #Or if a cluster is defined
#' cl <- makeCluster(12)
#' clusterCall(cl, function() library(tidyverse))
#' clusterCall(cl, function() library(optimr))
#' clusterCall(cl, function() library(SkeweDF))
#' parameters <- local_fit_function(param_bounds = param_bounds, data = obs_data, model_fn = Kolmogorov_Waring, clust = cl)
#' @export
local_fit_function <- function(param_bounds, data, model_fn, weighted_rt = F, par_chunk = 100, par_chunk_size = 10, n_cores = 1, clust){
  
  data <- data / sum(data)
  
  standalone <- F
  
  param_bounds <<- param_bounds
  weighted_rt <<- weighted_rt
  model_fn <<- model_fn
  par_chunk_size <<- par_chunk_size
  
  if(missing(clust)){
    standalone <- T
    clust <- makeCluster(n_cores)
    clusterExport(clust, varlist = c('data','param_bounds','weighted_rt', 'par_chunk_size'))
    clusterCall(clust, function() library(tidyverse))
    clusterCall(clust, function() library(optimr))
    clusterCall(clust, function() library(SkeweDF))
  }

  
  clusterExport(clust, varlist = c('param_bounds','weighted_rt', 'par_chunk_size','model_fn'))
  
  parameters <- parLapply(clust,1:par_chunk, function(q){
    par_mat <- lapply(1:par_chunk_size, function(i){
      par_vec <- c(rep(0, length(param_bounds)))
      for(i in 1:length(param_bounds)){
        par_vec[i] <- sample(param_bounds[[i]], size = 1, replace = T)
      }
      return(par_vec)
    }) %>% as.data.frame() %>% t() %>% as.matrix()
    
    param_lower <- c(rep(0, length(param_bounds)))
    param_upper <- c(rep(0, length(param_bounds)))
    
    for(i in 1:length(param_bounds)){
      param_lower[i] <- min(param_bounds[[i]])
      param_upper[i] <- max(param_bounds[[i]])
    }
    
    fn_parameters <- multistart(par = par_mat, fn = psi_criterion_function,method = 'L-BFGS-B', #control = list(fnscale = -1),
                                lower = param_lower, upper = param_upper, 
                                data = data, model_fn = model_fn, pmf_weight = 0.0, weighted_rt = weighted_rt); 
    
    colnames(fn_parameters)[1:(length(param_bounds)+1)] <- c(formalArgs(model_fn)[-1],'Psi_RTCDF')
    fn_parameters$Psi_RTCDF <- fn_parameters$Psi_RTCDF * -1
    return(fn_parameters)
  }) %>% bind_rows()

  print('Complete')
  
  if(standalone){
    stopCluster(clust)
  }
  
  parameters <- arrange(parameters, desc(Psi_RTCDF))
  return(parameters)

}

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
#' @examples
#' obs_data <- c(16362,7910,4563,2880,1992,1242,936,734,553,346,253,296,185,210,120,121,119,94,87,55,47,60,67,26,24,54,32,31,16,19,19,20,11,22,21,25,8,8,28,9,4,4,2,8,6,5,1,9,0,6,2,4,2,1,0,4,3,2,6,5,3,1,1,0,3,4,0,1,0,0,0,1,1,0,0,0,0,0,0,3,0,0,0,0,1,0,0,2,1,5,0,0,0,1,1,1,0,0,0,0,0,0,0,5,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4)
#' param_bounds <- list(seq(1,5,0.1),seq(1,5,0.1),seq(0.1,3,0.1),seq(0.1,3,0.1))
#' parameters <- local_fit_RGHD_ratio(param_bounds, obs_data, n_cores = 12)
#' 
#' #Or if a cluster is defined
#' cl <- makeCluster(12)
#' clusterCall(cl, function() library(tidyverse))
#' clusterCall(cl, function() library(optimr))
#' clusterCall(cl, function() library(SkeweDF))
#' parameters <- local_fit_RGHD_ratio(param_bounds, obs_data, clust = cl)
#' @export
local_fit_RGHD_ratio <- function(param_bounds, data, weighted_rt = F, par_chunk = 100, par_chunk_size = 10, n_cores = 1, clust){
  
  if(length(param_bounds) %% 2 != 0){
    stop('Must have even number of parameters!')
  }
  
  data <- data / sum(data)
  
  standalone <- F
  
  m <<- length(param_bounds) / 2
  param_bounds <<- param_bounds
  weighted_rt <<- weighted_rt
  
  if(missing(clust)){
    standalone <- T
    clust <- makeCluster(n_cores)
    clusterExport(clust, varlist = c('data','param_bounds','weighted_rt', 'par_chunk_size','m'))
    clusterCall(clust, function() library(tidyverse))
    clusterCall(clust, function() library(optimr))
    clusterCall(clust, function() library(SkeweDF))
  }
  
  
  clusterExport(clust, varlist = c('param_bounds','weighted_rt', 'par_chunk_size','m'))
  
  parameters <- parLapply(clust,1:par_chunk, function(q){
    par_mat <- lapply(1:par_chunk_size, function(i){
      par_vec <- c(rep(0, length(param_bounds)))
      for(i in 1:length(param_bounds)){
        par_vec[i] <- sample(param_bounds[[i]], size = 1, replace = T)
      }
      return(par_vec)
    }) %>% as.data.frame() %>% t() %>% as.matrix()
    
    param_lower <- c(rep(0, length(param_bounds)))
    param_upper <- c(rep(0, length(param_bounds)))
    
    for(i in 1:length(param_bounds)){
      param_lower[i] <- min(param_bounds[[i]])
      param_upper[i] <- max(param_bounds[[i]])
    }
    
    fn_parameters <- multistart(par = par_mat,fn = psi_criterion_RGHD_ratio,method = 'L-BFGS-B',
                                lower = param_lower, upper = param_upper, 
                                data = data, m = m, pmf_weight = 0.0, weighted_rt = weighted_rt);
    
    for(i in 1:m){
      colnames(fn_parameters)[i] <- paste0('r',i)
      colnames(fn_parameters)[i+m] <- paste0('q',i)
    }
    colnames(fn_parameters)[m+m+1] <- 'Psi_RTCDF'
    fn_parameters$Psi_RTCDF <- fn_parameters$Psi_RTCDF * -1
    return(fn_parameters)
  }) %>% bind_rows()
  
  print('Complete')
  
  if(standalone){
    stopCluster(clust)
  }
  
  parameters <- arrange(parameters, desc(Psi_RTCDF))
  return(parameters)
  
}

#' Global optimization of a given function given empirical data and parameter bounds
#'
#' This function generates a single set of optimized parameters and Psi Criterion for a given function within specified starting parameter bounds. This function uses a modified grid search method for optimization
#' @param param_bounds A list of sequences which indicate space where parameters should be generated and fit
#' @param data Vector of observed values 
#' @param model_fn Function of theoretical model to be used. For example, for Generalized_Yule(n, rho, alpha), model_fn <- Generalied_Yule
#' @param iter Integer indicating number of iterations to run grid search. Increasing iterations will increase decimal point precision of output parameters.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param n_cores Integer used to indicate number of cores to be used for this function if a socket cluster object is not defined.
#' @param clust socket cluster object from 'parallel::makeCluster()'. This is used if you have already generated a socket cluster object and would like to run this functoin on it. If no object is defined, one will be made for this function call.
#' @examples
#' obs_data <- c(16362,7910,4563,2880,1992,1242,936,734,553,346,253,296,185,210,120,121,119,94,87,55,47,60,67,26,24,54,32,31,16,19,19,20,11,22,21,25,8,8,28,9,4,4,2,8,6,5,1,9,0,6,2,4,2,1,0,4,3,2,6,5,3,1,1,0,3,4,0,1,0,0,0,1,1,0,0,0,0,0,0,3,0,0,0,0,1,0,0,2,1,5,0,0,0,1,1,1,0,0,0,0,0,0,0,5,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4)
#' param_bounds <- list(seq(0.1,5,0.1),seq(0.1,5,0.1),seq(0.01,1,0.01))
#' parameters <- global_fit_function(param_bounds = param_bounds, data = obs_data, model_fn = Kolmogorov_Waring, iter = 3, n_cores = 12)
#' 
#' #Or if a cluster is defined
#' cl <- makeCluster(12)
#' clusterCall(cl, function() library(tidyverse))
#' clusterCall(cl, function() library(optimr))
#' clusterCall(cl, function() library(SkeweDF))
#' parameters <- global_fit_function(param_bounds = param_bounds, data = obs_data, model_fn = Kolmogorov_Waring, iter = 3, clust = cl)
#' @export
global_fit_function <- function(param_bounds, data, model_fn, iter = 1, weighted_rt = F, n_cores = 1, clust){
  
  standalone <- F
  
  if(missing(clust)){
   standalone <- T
   clust <- makeCluster(n_cores)
   clusterExport(clust, varlist = c('data'))
   clusterCall(clust, function() library(tidyverse))
   clusterCall(clust, function() library(optimr))
   clusterCall(clust, function() library(SkeweDF))
  }
  
  n_parameters <<- length(param_bounds)
  
  par_mat <- expand.grid(param_bounds)
  par_mat$n <- length(data)
  par_mat <- par_mat[c(length(par_mat),1:(length(par_mat)-1))]
  
  output <<- par_mat
  
  data <- data / sum(data)
  if(weighted_rt){
    right_cdf <<- weighted_right_tail_cdf(data)
  }else{
    right_cdf <<- right_tail_cdf(data)
  }
  
  print(paste('Parameter space generated - # parameters:', n_parameters))
  
  clusterExport(clust, varlist = c('data','output','n_parameters','right_cdf'))
  clusterCall(clust, function() model_function <<- model_fn)
  
  
  
  criterion_list <- parLapply(clust, 1:nrow(output),function(i){ 
    model <- invoke(model_function, output[i,] %>% unlist() %>% unname())
    if(weighted_rt){
      model_right_cdf <- weighted_right_tail_cdf(model)
    }
    else{
      model_right_cdf <- right_tail_cdf(model)
    }
    return(psi_criterion(right_cdf, model_right_cdf, n_parameters))
  }) %>% unlist()

  print('Iteration: 1')
  
  
  output$criterion <- criterion_list
  output <- output %>% filter(!is.na(criterion))
  
  if(iter >= 2){
    for(q in 1:(iter-1)){
      bVars <- output[output$criterion == max(output$criterion),][1,] %>% unlist()
      bVars <- bVars[2:(length(bVars)-1)]
      
      bVars_decimal <- gsub("^.*\\.","",  bVars %>% as.character()) %>% nchar()
      
      bVars_decimal <- bVars_decimal + 1
      
      seq_list <- 1:length(bVars) %>% as.list()
      
      for(i in 1:length(seq_list)){
        seq_list[[i]] <- seq(bVars[i] - 5 * (10 ^ -bVars_decimal[i]),bVars[i] + 5 * (10 ^ -bVars_decimal[i]),by = 10 ^ -bVars_decimal[i]) %>% round(bVars_decimal[i])
        
      }
      
      output <- expand.grid(seq_list)
      output$n <- length(data)
      output <- output[,c(length(output),1:(length(output)-1))]
      colnames(output) <- formalArgs(model_fn)
      
      clusterExport(clust, varlist = c('output'))
      criterion_list <- parLapply(clust, 1:nrow(output),function(i){ 
        model <- invoke(model_function, output[i,] %>% unlist() %>% unname())
        if(weighted_rt){
          model_right_cdf <- weighted_right_tail_cdf(model)
        }
        else{
          model_right_cdf <- right_tail_cdf(model)
        }
        return(psi_criterion(right_cdf, model_right_cdf, n_parameters))
      }) %>% unlist()
      
      output$criterion <- criterion_list
      output <- output %>% filter(!is.na(criterion))
      
      print(paste('Iteration:', q+1))
    }
  }
  
  print('Complete')
  
  if(standalone){
    stopCluster(clust)
  }
  
  params <- output[output$criterion == max(output$criterion),][1,]
  return(params)
}

#' Global optimization of the 2m-RGHD function given empirical data, r bounds, and q/r bounds.
#'
#' This function generates a single set of optimized parameters and Psi Criterion for a given function within specified starting parameter bounds. This function uses Limited Memory BFGS as it's gradient descent algorithm.
#' @param param_bounds A list of sequences which indicate space where parameters should be generated and fit
#' @param data Vector of observed values 
#' @param iter Integer indicating number of iterations to run grid search. Increasing iterations will increase decimal point precision of output parameters.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param n_cores Integer used to indicate number of cores to be used for this function if a socket cluster object is not defined.
#' @param clust socket cluster object from 'parallel::makeCluster()'. This is used if you have already generated a socket cluster object and would like to run this functoin on it. If no object is defined, one will be made for this function call.
#' @examples
#' obs_data <- c(16362,7910,4563,2880,1992,1242,936,734,553,346,253,296,185,210,120,121,119,94,87,55,47,60,67,26,24,54,32,31,16,19,19,20,11,22,21,25,8,8,28,9,4,4,2,8,6,5,1,9,0,6,2,4,2,1,0,4,3,2,6,5,3,1,1,0,3,4,0,1,0,0,0,1,1,0,0,0,0,0,0,3,0,0,0,0,1,0,0,2,1,5,0,0,0,1,1,1,0,0,0,0,0,0,0,5,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4)
#' param_bounds <- list(seq(1,5,0.1),seq(1,5,0.1),seq(0.1,3,0.1),seq(0.1,3,0.1))
#' parameters <- global_fit_RGHD_ratio(param_bounds, obs_data, iter = 3, n_cores = 12)
#' 
#' #Or if a cluster is defined
#' cl <- makeCluster(12)
#' clusterCall(cl, function() library(tidyverse))
#' clusterCall(cl, function() library(optimr))
#' clusterCall(cl, function() library(SkeweDF))
#' parameters <- global_fit_RGHD_ratio(param_bounds, obs_data, iter = 3, clust = cl)
#' @export
global_fit_RGHD_ratio <- function(param_bounds, data, iter, weighted_rt = F, n_cores = 1, clust){
  
  if(length(param_bounds) %% 2 != 0){
    stop('Must have even number of parameters!')
  }
  
  standalone <- F
  
  if(missing(clust)){
    standalone <- T
    clust <- makeCluster(n_cores)
    clusterExport(clust, varlist = c('data'))
    clusterCall(clust, function() library(tidyverse))
    clusterCall(clust, function() library(optimr))
    clusterCall(clust, function() library(SkeweDF))
  }
  
  m <- length(param_bounds) / 2
  n_parameters <<- m * 2
  
  par_mat <- expand.grid(param_bounds)
  par_mat$n <- length(data)
  par_mat <- par_mat[c(length(par_mat),1:(length(par_mat)-1))]
  
  output <<- par_mat
  
  data <- data / sum(data)
  if(weighted_rt){
    right_cdf <<- weighted_right_tail_cdf(data)
  }else{
    right_cdf <<- right_tail_cdf(data)
  }
  
  print(paste('Parameter space generated - # parameters:', n_parameters))
  
  clusterExport(clust, varlist = c('data','output','n_parameters','right_cdf'))
  
  criterion_list <<- parLapply(clust, 1:nrow(output),function(i){ 
    model <- RGHD(length(data), m, unlist(output[i,1:m]), unlist(output[i,1:m]) * unlist(output[i,(m+1):(m+m)]))
    if(weighted_rt){
      model_right_cdf <- weighted_right_tail_cdf(model)
    }
    else{
      model_right_cdf <- right_tail_cdf(model)
    }
    return(psi_criterion(right_cdf, model_right_cdf, n_parameters))
  }) %>% unlist()
  
  print('Iteration: 1')
  
  
  output$criterion <<- criterion_list
  output <<- output %>% filter(!is.na(criterion))
  
  if(iter >= 2){
    for(q in 1:(iter-1)){
      bVars <- output[output$criterion == max(output$criterion),][1,] %>% unlist()
      bVars <- bVars[2:(length(bVars)-1)]
      
      bVars_decimal <- gsub("^.*\\.","",  bVars %>% as.character()) %>% nchar()
      
      bVars_decimal <- bVars_decimal + 1
      
      seq_list <- 1:length(bVars) %>% as.list()
      
      for(i in 1:length(seq_list)){
        seq_list[[i]] <- seq(bVars[i] - 5 * (10 ^ -bVars_decimal[i]),bVars[i] + 5 * (10 ^ -bVars_decimal[i]),by = 10 ^ -bVars_decimal[i]) %>% round(bVars_decimal[i])
        
      }
      
      output <<- expand.grid(seq_list)
      output$n <<- length(data)
      output <<- output[,c(length(output),1:(length(output)-1))]
      
      clusterExport(clust, varlist = c('output'))
      criterion_list <- parLapply(clust, 1:nrow(output),function(i){ 
        model <- RGHD(length(data), m, unlist(output[i,1:m]), unlist(output[i,1:m]) * unlist(output[i,(m+1):(m+m)]))
        if(weighted_rt){
          model_right_cdf <- weighted_right_tail_cdf(model)
        }
        else{
          model_right_cdf <- right_tail_cdf(model)
        }
        return(psi_criterion(right_cdf, model_right_cdf, n_parameters))
      }) %>% unlist()
      
      output$criterion <- criterion_list
      output <- output %>% filter(!is.na(criterion))
      
      print(paste('Iteration:', q+1))
    }
  }
  
  print('Complete')
  
  if(standalone){
    stopCluster(clust)
  }
  
  for(i in 1:m){
    colnames(output)[i+1] <- paste0('r',i)
    colnames(output)[i+1+m] <- paste0('q',i)
  }
  
  params <- output[output$criterion == max(output$criterion),][1,]
  return(params)
}

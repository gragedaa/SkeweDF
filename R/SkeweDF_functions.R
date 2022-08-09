## usethis namespace: start
#' @useDynLib SkeweDF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats sd qt lm cor ks.test
#' @importFrom zipfR Ibeta
#' @importFrom dplyr %>% filter desc arrange bind_rows
#' @importFrom stringr str_replace
#' @importFrom parallel makeCluster clusterExport clusterCall stopCluster parLapply
#' @importFrom purrr invoke
#' @importFrom optimr multistart
#' @importFrom methods formalArgs
#' @importFrom matrixStats colSds colMedians
#' @importFrom grDevices dev.off png rgb
#' @importFrom graphics abline plot points text
#' @importFrom utils write.table
## usethis namespace: end
NULL


#' Get Mean Confidence Interval Function
#'
#' This function generates a vector of confidence interval based on mean of data.
#' @param data Data to get confidence interval from
#' @param alpha Alpha for confidence interval calculation
#' @export
get_CI <- function(data, alpha){
  std <- sd(data)
  mean <- mean(data)
  sample_size <- length(data)
  mean_error <- qt(1-alpha, sample_size-1) * std / sqrt(sample_size)

  return(c(mean - mean_error, mean + mean_error))
}

#' Get Median Confidence Interval Function
#'
#' This function generates a vector of ranked 95% confidence interval based on median of data.
#' @param data Data to get confidence interval from
#' @export
get_median_CI <- function(data){
  n <- length(data)
  lower_rank <- (n/2) - (1.96 * sqrt(n) / 2)
  upper_rank <- 1 + (n/2) + (1.96 * sqrt(n) / 2)
  lower_rank <- lower_rank %>% round(0)
  upper_rank <- upper_rank %>% round(0)

  data <- sort(data)
  return(c(data[lower_rank], data[upper_rank]))
}

#' Exponential Distribution Function
#'
#' This function generates a vector of n length of the Exponential distribution with parameters a and b.
#' @param n Length of vector to be generated.
#' @param a Parameter of the Exponential distribution function
#' @param b Parameter of the Exponential distribution function
#' @examples
#' Exponential(100, 10000, 0.8)
#' @export
Exponential <- function(n, a , b){
  p <- a * exp(-b * (1:n))
  return(p)
}

#' Lorentzian Distribution Function calculation
#'
#' This function calculates value of Lorentzian function at x
#' @param x Index of function
#' @param gamma Parameter of the Lorenzian distribution function
#' @param x0 Parameter of the Lorenzian distribution function indicating center of function
#' @param c Parameter of the Lorenzian distribution function indicating center of function
#' @examples
#' Lorentzian_calc(5, 5.5, 6, 2)
#' @export
Lorentzian_calc <- function(x, gamma, x0, c){

  a <- abs(x- x0) ^ c

  b <- gamma / 2

  output <- b / (a + b^c)
  return(output)
}

#' Lorentzian Distribution Function
#'
#' This function generates a vector of n length of the Lorentzian distribution
#' @param n Length of vector to be generated.
#' @param gamma Parameter of the Lorenzian distribution function
#' @param x0 Parameter of the Lorenzian distribution function indicating center of function
#' @param c Parameter of the Lorenzian distribution function indicating center of function
#' @examples
#' Lorentzian_calc(5, 5.5, 6, 2)
#' @export
Lorentzian <- function(n, gamma, x0, c){
  out <- lapply(1:n, function(i){
         return(Lorentzian_calc(i, gamma, x0, c))
     })
  out <- unlist(out)
  out <- out / sum(out)
  return(out)
}

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

#' Weighted Left-Tail Cumulative Distribution Function
#'
#' This function generates a vector of the weighted left-tail cumulative distribution function of a given vector of values. The weight of of each variable is determined by its position in the vector. For example, with a vector of length 5, element 1 will have weight 5/(5+4+3+2+1). Element 1 will have weight 5/(5+4+3+2+1)
#' @param x Length of vector to be generated.
#' @examples
#' x <- c(1,2,3,4,5)
#' weighted_left_tail_cdf(x)
#' @export
weighted_left_tail_cdf <- function(x){
  output_cdf <- x * length(x):1
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

#' Local optimization of a given function given empirical data and parameter bounds
#'
#' This function generates a table of optimized parameters and Psi Criterion for a given function within specified starting parameter bounds. This function uses Limited Memory BFGS as it's gradient descent algorithm.
#' @param param_bounds A list of sequences which indicate space where parameters should be generated and fit
#' @param data Vector of observed values
#' @param model_fn_name Character vector indicating name of function of theoretical model to be used. For example, for Generalized_Yule(n, rho, alpha), model_fn_name <- 'Generalized_Yule'
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param par_chunk Integer used to indicate number of optimization chunks to be run. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param par_chunk_size Integer used to indicate number of starting parameters to be generated and optimized in a given chunk. Total number of rows in the output table = par_chunk * par_chunk_size
#' @param n_cores Integer used to indicate number of cores to be used for this function if a socket cluster object is not defined.
#' @param clust socket cluster object from 'parallel::makeCluster()'. This is used if you have already generated a socket cluster object and would like to run this functoin on it. If no object is defined, one will be made for this function call.
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @export
local_fit_function <- function(param_bounds, data, model_fn_name, weighted_rt = FALSE, pmf_weight = 0.0, par_chunk = 100, par_chunk_size = 10, n_cores = 1, clust, left_trunc = 1, right_trunc = left_trunc+length(data)-1){

  skeweDF_gVar_data <- data / sum(data)

  standalone <- FALSE

  skeweDF_gVar_param_bounds <- param_bounds
  skeweDF_gVar_weighted_rt <- weighted_rt
  skeweDF_gVar_model_fn <- get(model_fn_name)
  skeweDF_gVar_par_chunk_size <- par_chunk_size
  skeweDF_gVar_left_trunc <- left_trunc
  skeweDF_gVar_right_trunc <- right_trunc

  if(missing(clust)){
    standalone <- TRUE
    clust <- makeCluster(n_cores)
    clusterExport(clust, varlist = c('skeweDF_gVar_data','skeweDF_gVar_param_bounds','skeweDF_gVar_weighted_rt', 'skeweDF_gVar_par_chunk_size','skeweDF_gVar_model_fn','skeweDF_gVar_left_trunc','skeweDF_gVar_right_trunc', 'pmf_weight'), envir = environment())
    clusterCall(clust, function() library(SkeweDF))
  }


  clusterExport(clust, varlist = c('skeweDF_gVar_data','skeweDF_gVar_param_bounds','skeweDF_gVar_weighted_rt', 'skeweDF_gVar_par_chunk_size','skeweDF_gVar_model_fn','skeweDF_gVar_left_trunc','skeweDF_gVar_right_trunc', 'pmf_weight'), envir = environment())

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

    fn_parameters <- multistart(parmat = par_mat, fn = psi_criterion_function,method = 'L-BFGS-B', #control = list(fnscale = -1),
                                lower = param_lower, upper = param_upper,
                                data = skeweDF_gVar_data, model_fn = skeweDF_gVar_model_fn, pmf_weight = pmf_weight, weighted_rt = skeweDF_gVar_weighted_rt,
                                left_trunc = skeweDF_gVar_left_trunc, right_trunc = skeweDF_gVar_right_trunc);

    colnames(fn_parameters)[1:(length(skeweDF_gVar_param_bounds)+1)] <- c(formalArgs(skeweDF_gVar_model_fn)[-1],'Psi_RTCDF')
    fn_parameters$Psi_RTCDF <- fn_parameters$Psi_RTCDF * -1
    return(fn_parameters)
  }) %>% bind_rows()

  print('Complete')

  if(standalone){
    stopCluster(clust)
  }

  parameters <- arrange(parameters, desc(parameters$Psi_RTCDF))
  parameters$fevals <- NULL
  parameters$gevals <- NULL
  parameters$convergence <- NULL

  parameters <- parameter_post_processing(parameters, model_fn_name, skeweDF_gVar_data)


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


#' Parameter Optimization Helper Function
#'
#' This function adds in additional columns to the optimized parameter output dataframe
#' @param parameter_df Output dataframe of optimized parameters using local algorithm
#' @param model_fn_name Character vector used to indicate name of model function used for optimization
#' @param data Vector of observed values
#' @export
parameter_post_processing <- function(parameter_df, model_fn_name, data){

  #Psi PDF, Psi weighted RTCDF, Error sum of squares pdf, Esos cdf, esos wtrtcdf, r2, ks pdf, ks rtcdf, ks weighted rtcdf

  n_parameters <- length(parameter_df) - 1

  model_fn <- get(model_fn_name)
  model_list <- vector('list',nrow(parameter_df))

  if(!(model_fn_name == 'RGHD')){

    model_list <- lapply(1:length(model_list), function(i){
      output <- invoke(model_fn, c(length(data), parameter_df[i,1:n_parameters]) %>% unlist() %>% unname())
      output <- output / sum(output)
      return(output)
    })

  }
  else{
    m <- n_parameters / 2

    model_list <- lapply(1:length(model_list), function(i){
      output <- RGHD(length(data), m, c(parameter_df[i,1:m] %>% unlist() %>% unname()), c(parameter_df[i,(m+1):(m+m)] %>% unlist() %>% unname()))
      output <- output / sum(output)
      return(output)
    })
  }

  if(model_fn_name == 'Kolmogorov_Waring'){
    #get p0
    parameter_df$p0 <- 0
    for(i in 1:nrow(parameter_df)){
      parameter_df$p0[i] <- get_p0(parameter_df[i,-c(length(parameter_df)-1, length(parameter_df))] %>% unlist(),model_fn_name)
    }

    #calculate a/b ratio
    parameter_df$ab_ratio <- parameter_df$a / parameter_df$b
  }
  else if(model_fn_name == 'RGHD'){
    m <- n_parameters / 2

    #get p0
    parameter_df$p0 <- 0
    for(i in 1:nrow(parameter_df)){
      parameter_df$p0[i] <- get_p0(parameter_df[i,-c(length(parameter_df)-1, length(parameter_df))] %>% unlist(),'RGHD')
    }

    #get ratios
    for(i in 1:m){
      parameter_df[paste0('r',i,'q',i,'_ratio')] <-  parameter_df[i] / parameter_df[i+m]
    }

      if( m > 1){
        #reorder such that r1q1_ratio > r2q2_ratio > r3q3_ratio ...

        for(i in 1:nrow(parameter_df)){
          tmp_params <- parameter_df[i,1:(m+m)] %>% unlist()
          tmp_ratios <- parameter_df[i,(length(parameter_df)-m+1):(length(parameter_df))] %>% unlist()
          ranked_ratios <- rank(-tmp_ratios)

          output_params <- tmp_params
          output_ratios <- tmp_ratios

          for( r in 1:length(ranked_ratios)){
            ratio_rank <- ranked_ratios[r]
            output_params[c(ratio_rank, ratio_rank + m)] <- tmp_params[c(r,r+m)]
            output_ratios[ratio_rank] <- tmp_ratios[r]
          }

          parameter_df[i,1:(m+m)] <- output_params
          parameter_df[i,(length(parameter_df)-m+1):(length(parameter_df))] <- output_ratios
      }

    }
  }

  parameter_df$Psi_PDF <- lapply(1:length(model_list), function(i){
    psi_criterion(data / sum(data), model = model_list[[i]], n_parameters = n_parameters)
  }) %>% unlist()

  parameter_df$Psi_weighted_RTCDF <- lapply(1:length(model_list), function(i){
    psi_criterion(weighted_right_tail_cdf(data / sum(data)), model = weighted_right_tail_cdf(model_list[[i]]), n_parameters = n_parameters)
  }) %>% unlist()

  parameter_df$Error_sum_of_squares_PDF <- lapply(1:length(model_list), function(i){
    one <- data / sum(data)
    two <- model_list[[i]]

    output <-  sum((one - two) ^ 2)

    return(output)
  }) %>% unlist()

  parameter_df$Error_sum_of_squares_RTCDF <- lapply(1:length(model_list), function(i){
    one <- data / sum(data)
    two <- model_list[[i]]

    output <-  sum((right_tail_cdf(one) - right_tail_cdf(two)) ^ 2)

    return(output)
  }) %>% unlist()

  parameter_df$Error_sum_of_squares_weighted_RTCDF <- lapply(1:length(model_list), function(i){
    one <- data / sum(data)
    two <- model_list[[i]]

    output <-  sum((weighted_right_tail_cdf(one) - weighted_right_tail_cdf(two)) ^ 2)

    return(output)
  }) %>% unlist()

  parameter_df$R_squared <- lapply(1:length(model_list), function(i){
    output <-  cor(data, model_list[[i]]) ^ 2

    return(output)
  }) %>% unlist()

  parameter_df$KS_PDF <- lapply(1:length(model_list), function(i){
    output <-  ks.test(data / sum(data), model_list[[i]])$statistic

    return(output)
  }) %>% unlist()

  parameter_df$KS_RTCDF <- lapply(1:length(model_list), function(i){
    output <-  ks.test(right_tail_cdf(data / sum(data)), right_tail_cdf(model_list[[i]]))$statistic

    return(output)
  }) %>% unlist()

  parameter_df$KS_weighted_RTCDF <- lapply(1:length(model_list), function(i){
    output <-  ks.test(weighted_right_tail_cdf(data / sum(data)), weighted_right_tail_cdf(model_list[[i]]))$statistic

    return(output)
  }) %>% unlist()


  return(parameter_df)
}

#' Label Coordinate Calculate Helper Function
#'
#' This function calculates coordinates for a plot given x and y bounds and location represented as percentage of plot area
#' @param x_lower_bound Numeric lowest value of x axis
#' @param x_upper_bound Numeric highest value of x axis
#' @param y_lower_bound Numeric lowest value of y axis
#' @param y_upper_bound Numeric highest value of y axis
#' @param x_buffer Numeric indicating location on x axis (0 - 1)
#' @param y_buffer Numeric indicating location on y axis (0 - 1)
#' @param log_scale_x Boolean indicating if x axis is log scale
#' @param log_scale_y Boolean indicating if y axis is log scale
#' @export
calculate_label_coords <- function(x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound, x_buffer = 0.5, y_buffer = 0.5, log_scale_x = FALSE, log_scale_y = FALSE){
  output_x <- 0
  output_y <- 0

  if(log_scale_x){
    x_lower_bound <- log(x_lower_bound,10)
    x_upper_bound <- log(x_upper_bound,10)
  }
  if(log_scale_y){
    y_lower_bound <- log(y_lower_bound,10)
    y_upper_bound <- log(y_upper_bound,10)
  }

  diff_x <- x_upper_bound - x_lower_bound
  diff_y <- y_upper_bound - y_lower_bound

  output_x <- x_lower_bound + (x_buffer * diff_x)
  output_y <- y_lower_bound + (y_buffer * diff_y)

  if(log_scale_x){
    output_x <- 10 ^ output_x
  }
  if(log_scale_y){
    output_y <- 10 ^ output_y
  }

  return(c(output_x, output_y))
}

#' Plot Model Helper Function
#'
#' This function generates various plots of empirical data and models
#' @param title Character vector indicating title of the empirical dataset, this will be present on every plot, this also determines the name of the folder where plots will be
#' @param model_fn_name Character vector used to indicate name of model function used for optimization
#' @param data Vector of observed values
#' @param parameter_df Data frame of optimized parameters and other model function values (p0, Psi, etc)
#' @param n_parameters Int of number of parameters used in model funciton
#' @param plot_folder_name Character vector indicating folder or directory name to be used when outputting plot images
#' @param xlab Character vector indicating x axis label of plots, indicates what the random variable is
#' @param left_trunc Int indicating starting index of model function used for optimization
#' @export
plot_model <- function(title, model_fn_name, data, parameter_df, n_parameters, plot_folder_name, xlab, left_trunc = 1){

  all_parameter_df <- parameter_df
  if(model_fn_name == 'Kolmogorov_Waring' | model_fn_name == 'RGHD'){
    p0_med <- get_median_CI(parameter_df$p0)
    parameter_df <- parameter_df[parameter_df$p0 >= p0_med[1] & parameter_df$p0 <= p0_med[2],]
  }

  bootstrap_n <- nrow(parameter_df)

  model_fn <- get(model_fn_name)
  model_list <- vector('list',nrow(parameter_df))

  if(!(model_fn_name == 'RGHD')){
    for(i in 1:length(model_list)){
      model_list[[i]] <- invoke(model_fn, c(length(data), parameter_df[i,1:n_parameters]) %>% unlist() %>% unname())
      model_list[[i]] <- model_list[[i]] /sum(model_list[[i]][left_trunc:length(data)])
    }
  }
  else{
    m <- n_parameters / 2
    for(i in 1:length(model_list)){
      model_list[[i]] <- RGHD(length(data), m, c(parameter_df[i,1:m] %>% unlist() %>% unname()), c(parameter_df[i,(m+1):(m+m)] %>% unlist() %>% unname()))
      model_list[[i]] <- model_list[[i]] /sum(model_list[[i]][left_trunc:length(data)])
    }
  }



  model <- model_list[[1]]


  #replace parameter df with tmp
  if(!(model_fn_name == 'RGHD')){
    fn_name <- str_replace(model_fn_name,'_',' ')

    plot_label <- paste0(names(parameter_df[1]),': ', signif(parameter_df[1,1], digits=3))
    if(n_parameters > 1){
      for(i in 2:n_parameters){
        plot_label <- paste0(plot_label,'\n',names(parameter_df[i]),': ', signif(parameter_df[1,i], digits=3))
      }
    }



  }
  else{
    fn_name <- paste0('2m-RGHD (m=', n_parameters/2, ')')

    plot_label <- paste0(names(parameter_df[1]),': ', signif(parameter_df[1,1], digits=3),
                         ' ',names(parameter_df[(n_parameters/2)+1]),': ', signif(parameter_df[1,(n_parameters/2)+1], digits=3))
    if(n_parameters > 2){
      for(i in 2:(n_parameters/2)){
        plot_label <- paste0(plot_label,'\n', names(parameter_df[i]),': ', signif(parameter_df[1,i], digits=3),
                             ' ',names(parameter_df[(n_parameters/2)+i]),': ', signif(parameter_df[1,(n_parameters/2)+i], digits=3))
      }
    }

  }

  if(model_fn_name == 'Kolmogorov_Waring' | model_fn_name == 'RGHD'){
    plot_label <- paste0(plot_label,'\np0: ', signif(parameter_df[1,'p0'], digits=5))
  }

  plot_label <- paste0(plot_label,'\nPsi_RTCDF: ', signif(parameter_df[1,'Psi_RTCDF'], digits=5))



  dir.create(paste0(plot_folder_name,'/',fn_name))


  plot_y_floor = 10 ^ (min(data[data != 0]) %>% log(10) %>% floor())

  tmp <- calculate_label_coords(1, length(data), min(data[data != 0]), max(data), x_buffer = 0, y_buffer = 0.25, log_scale_y = TRUE)
  pmf_text_x <- tmp[1]
  pmf_text_y <- tmp[2]

  tmp <- calculate_label_coords(1, length(data), 0, sum(data[left_trunc:length(data)]), x_buffer = 0.95, y_buffer = 0.85)
  cdf_text_x <- tmp[1]
  cdf_text_y <- tmp[2]


  png(paste0(plot_folder_name,'/',fn_name,'/000.png'), width = 2000, height = 2000, res = 300)
  plot(1:length(data), model * sum(data[left_trunc:length(data)]), log = 'xy',pch = 16, col = 'red', ylim = c(plot_y_floor,max(max(model), max(data))),
       main = paste0(title,'\n',fn_name),
       xlab = xlab,
       ylab = 'Frequency')
  points(1:length(data), data)
  text(pmf_text_x,pmf_text_y,pos = 4,labels = plot_label)
  dev.off()

  png(paste0(plot_folder_name,'/',fn_name,'/001.png'), width = 2000, height = 2000, res = 300)
  plot(1:length(data), model * sum(data[left_trunc:length(data)]),pch = 16, log = 'xy', col = 'red',
       main = paste0(title,'\n',fn_name),
       xlab = xlab,
       ylab = 'Frequency')
  points(1:length(data), data)
  text(pmf_text_x,pmf_text_y,pos = 4,labels = plot_label)
  dev.off()

  if(length(model_list) > 1){
    png(paste0(plot_folder_name,'/',fn_name,'/002.png'), width = 2000, height = 2000, res = 300)
    plot(1:length(data), model * sum(data[left_trunc:length(data)]), log = 'xy',pch = 16, col = 'red', ylim = c(plot_y_floor,max(max(model), max(data))),
         main = paste0(title,'\n',fn_name,'\nTop 5%'),
         xlab = xlab,
         ylab = 'Frequency')
    for(i in 1:(bootstrap_n * 0.05)){
      points(1:length(data), model_list[[i]] * sum(data[left_trunc:length(data)]), col = 'red',pch = 16)
    }
    points(1:length(data), data)
    points(1:length(data), model * sum(data[left_trunc:length(data)]), col = 'blue', pch = 16, cex = 0.5)
    text(pmf_text_x,pmf_text_y,pos = 4,labels = paste0('Psi_RTCDF: ', signif(mean(parameter_df['Psi_RTCDF'] %>% unlist()), digits=3), ' +- ', signif(sd(parameter_df['Psi_RTCDF'] %>% unlist()), digits=3)))
    dev.off()
  }

  if(length(model_list) > 1){
    png(paste0(plot_folder_name,'/',fn_name,'/003.png'), width = 2000, height = 2000, res = 300)
    plot(1:length(data), model * sum(data[left_trunc:length(data)]),pch = 16, log = 'xy', col = 'red',
         main = paste0(title,'\n',fn_name,'\nTop 5%'),
         xlab = xlab,
         ylab = 'Frequency')
    for(i in 1:(bootstrap_n * 0.05)){
      points(1:length(data), model_list[[i]] * sum(data[left_trunc:length(data)]), col = 'red',pch = 16)
    }
    points(1:length(data), data)
    points(1:length(data), model * sum(data[left_trunc:length(data)]), col = 'blue', pch = 16, cex = 0.5)
    text(pmf_text_x,pmf_text_y,pos = 4,labels = paste0('Psi_RTCDF: ', signif(mean(parameter_df['Psi_RTCDF'] %>% unlist()), digits=3), ' +- ', signif(sd(parameter_df['Psi_RTCDF'] %>% unlist()), digits=3)))
    dev.off()
  }

  lr <- data %>% log(10) %>% as.data.frame()
  colnames(lr) <- 'emp_data'
  lr$model <- (model * sum(data[left_trunc:length(data)])) %>% log(10)
  lr <- lr[lr$emp_data != -Inf,]
  lr_sum <- lm(formula = model ~ emp_data, data = lr) %>% summary()

  range <- abs(abs(min(min(lr$emp_data),min(lr$model))) - abs(max(max(lr$emp_data),max(lr$model))))

  png(paste0(plot_folder_name,'/',fn_name,'/004.png'), width = 2000, height = 2000, res = 300)
  plot(lr$emp_data, lr$model, xlim = c(min(min(lr$emp_data),min(lr$model)), max(max(lr$emp_data),max(lr$model))), ylim = c(min(min(lr$emp_data),min(lr$model)), max(max(lr$emp_data),max(lr$model))),
       main = paste0(title,'\n',fn_name,'\nLog-Q-Q Plot'),
       xlab = paste0(title, ' Data'),
       ylab = fn_name)
  abline(lm(model ~ emp_data, data = lr))
  abline(a = lr_sum$coefficients[1,'Estimate'] - lr_sum$coefficients[1,'Std. Error'], b = lr_sum$coefficients[2,'Estimate'] - lr_sum$coefficients[2,'Std. Error'], col = 'red', lty = 'dashed')
  abline(a = lr_sum$coefficients[1,'Estimate'] + lr_sum$coefficients[1,'Std. Error'], b = lr_sum$coefficients[2,'Estimate'] + lr_sum$coefficients[2,'Std. Error'],col = 'red', lty = 'dashed')
  text(min(min(lr$emp_data),min(lr$model)) ,
       min(min(lr$emp_data),min(lr$model)) + (range * 0.75), pos = 4,
       labels = paste0('Intercept: ',signif(lr_sum$coefficients[1,'Estimate'], digits=3),' +- ',signif(lr_sum$coefficients[1,'Std. Error'], digits=3),
                       '\nSlope: ',signif(lr_sum$coefficients[2,'Estimate'], digits=3),' +- ',signif(lr_sum$coefficients[2,'Std. Error'], digits=3)))
  dev.off()

  png(paste0(plot_folder_name,'/',fn_name,'/005.png'), width = 2000, height = 2000, res = 300)
  plot(1:length(data), (model * sum(data[left_trunc:length(data)])) %>% right_tail_cdf(), log = 'x', xlim = rev(range(1:length(data))), col = 'red', pch = 16,
       main = paste0(title,'\n',fn_name,'\nRight-Tail CDF'),
       xlab = xlab,
       ylab = 'Cumulative Frequency')
  points(1:length(data), data %>% right_tail_cdf())
  text(cdf_text_x,cdf_text_y,pos = 4,labels = plot_label)
  dev.off()

  if(length(model_list) > 1){
    png(paste0(plot_folder_name,'/',fn_name,'/006.png'), width = 2000, height = 2000, res = 300)
    plot(1:length(data), (model * sum(data[left_trunc:length(data)])) %>% right_tail_cdf(), log = 'x', xlim = rev(range(1:length(data))), col = 'red', pch = 16,
         main = paste0(title,'\n',fn_name,'\nRight-Tail CDF Top 5%'),
         xlab = xlab,
         ylab = 'Cumulative Frequency')
    for(i in 1:(bootstrap_n * 0.05)){
      points(1:length(data), (model_list[[i]] * sum(data[left_trunc:length(data)])) %>% right_tail_cdf(), col = 'red',pch = 16)
    }
    points(1:length(data), data %>% right_tail_cdf())
    points(1:length(data), (model * sum(data[left_trunc:length(data)])) %>% right_tail_cdf(), col = 'blue', pch = 16, cex = 0.5)
    text(cdf_text_x,cdf_text_y,pos = 4,labels = paste0('Psi_RTCDF: ', signif(mean(parameter_df['Psi_RTCDF'] %>% unlist()), digits=3), ' +- ', signif(sd(parameter_df['Psi_RTCDF'] %>% unlist()), digits=3)))
    dev.off()
  }

  if(model_fn_name == 'Kolmogorov_Waring'){

    png(paste0(plot_folder_name,'/',fn_name,'/007.png'), width = 2000, height = 2000, res = 300)
    plot(all_parameter_df$p0, all_parameter_df$Psi_RTCDF,pch = 16, col = rgb(0,0,0,alpha = 0.1),
         main = paste0(title,'\n',fn_name,'\np0 vs Psi_RTCDF'),
         xlab = 'p0',
         ylab = 'Psi_RTCDF')
    abline(v = get_CI(all_parameter_df$p0, 0.05)[1], lty = 1)
    abline(v = get_CI(all_parameter_df$p0, 0.05)[2], lty = 1)
    abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[1], lty = 1)
    abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[2], lty = 1)
    abline(v = get_median_CI(all_parameter_df$p0)[1], lty = 2)
    abline(v = get_median_CI(all_parameter_df$p0)[2], lty = 2)
    abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[1], lty = 2)
    abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[2], lty = 2)
    dev.off()

    png(paste0(plot_folder_name,'/',fn_name,'/008.png'), width = 2000, height = 2000, res = 300)
    plot(all_parameter_df$ab_ratio, all_parameter_df$Psi_RTCDF,pch = 16, col = rgb(0,0,0,alpha = 0.1),
         main = paste0(title,'\n',fn_name,'\na/b vs Psi_RTCDF'),
         xlab = 'a/b',
         ylab = 'Psi_RTCDF')
    abline(v = get_CI(all_parameter_df$ab_ratio, 0.05)[1])
    abline(v = get_CI(all_parameter_df$ab_ratio, 0.05)[2])
    abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[1])
    abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[2])
    abline(v = get_median_CI(all_parameter_df$ab_ratio)[1], lty = 2)
    abline(v = get_median_CI(all_parameter_df$ab_ratio)[2], lty = 2)
    abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[1], lty = 2)
    abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[2], lty = 2)
    dev.off()

    png(paste0(plot_folder_name,'/',fn_name,'/009.png'), width = 2000, height = 2000, res = 300)
    plot(all_parameter_df$ab_ratio, all_parameter_df$p0,pch = 16, col = rgb(0,0,0,alpha = 0.1),
         main = paste0(title,'\n',fn_name,'\na/b vs p0'),
         xlab = 'a/b',
         ylab = 'p0')
    abline(v = get_CI(all_parameter_df$ab_ratio, 0.05)[1])
    abline(v = get_CI(all_parameter_df$ab_ratio, 0.05)[2])
    abline(h = get_CI(all_parameter_df$p0, 0.05)[1])
    abline(h = get_CI(all_parameter_df$p0, 0.05)[2])
    abline(v = get_median_CI(all_parameter_df$ab_ratio)[1], lty = 2)
    abline(v = get_median_CI(all_parameter_df$ab_ratio)[2], lty = 2)
    abline(h = get_median_CI(all_parameter_df$p0)[1], lty = 2)
    abline(h = get_median_CI(all_parameter_df$p0)[2], lty = 2)
    dev.off()
  }else if(model_fn_name == 'RGHD'){

    if(n_parameters >= 2){
      png(paste0(plot_folder_name,'/',fn_name,'/007.png'), width = 2000, height = 2000, res = 300)
      plot(all_parameter_df$p0, all_parameter_df$Psi_RTCDF,pch = 16, col = rgb(0,0,0,alpha = 0.1),
           main = paste0(title,'\n',fn_name,'\np0 vs Psi_RTCDF'),
           xlab = 'p0',
           ylab = 'Psi_RTCDF')
      abline(v = get_CI(all_parameter_df$p0, 0.05)[1], lty = 1)
      abline(v = get_CI(all_parameter_df$p0, 0.05)[2], lty = 1)
      abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[1], lty = 1)
      abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[2], lty = 1)
      abline(v = get_median_CI(all_parameter_df$p0)[1], lty = 2)
      abline(v = get_median_CI(all_parameter_df$p0)[2], lty = 2)
      abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[1], lty = 2)
      abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[2], lty = 2)
      dev.off()

      png(paste0(plot_folder_name,'/',fn_name,'/008.png'), width = 2000, height = 2000, res = 300)
      plot(all_parameter_df$r1q1_ratio, all_parameter_df$Psi_RTCDF,pch = 16, col = rgb(0,0,0,alpha = 0.1),
           main = paste0(title,'\n',fn_name,'\nr1/q1 vs Psi_RTCDF'),
           xlab = 'r1/q1',
           ylab = 'Psi_RTCDF')
      abline(v = get_CI(all_parameter_df$r1q1_ratio, 0.05)[1])
      abline(v = get_CI(all_parameter_df$r1q1_ratio, 0.05)[2])
      abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[1])
      abline(h = get_CI(all_parameter_df$Psi_RTCDF, 0.05)[2])
      abline(v = get_median_CI(all_parameter_df$r1q1_ratio)[1], lty = 2)
      abline(v = get_median_CI(all_parameter_df$r1q1_ratio)[2], lty = 2)
      abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[1], lty = 2)
      abline(h = get_median_CI(all_parameter_df$Psi_RTCDF)[2], lty = 2)
      dev.off()

      png(paste0(plot_folder_name,'/',fn_name,'/009.png'), width = 2000, height = 2000, res = 300)
      plot(all_parameter_df$r1q1_ratio, all_parameter_df$p0,pch = 16, col = rgb(0,0,0,alpha = 0.1),
           main = paste0(title,'\n',fn_name,'\nr1q1 vs p0'),
           xlab = 'r1q1',
           ylab = 'p0')
      abline(v = get_CI(all_parameter_df$r1q1_ratio, 0.05)[1])
      abline(v = get_CI(all_parameter_df$r1q1_ratio, 0.05)[2])
      abline(h = get_CI(all_parameter_df$p0, 0.05)[1])
      abline(h = get_CI(all_parameter_df$p0, 0.05)[2])
      abline(v = get_median_CI(all_parameter_df$r1q1_ratio)[1], lty = 2)
      abline(v = get_median_CI(all_parameter_df$r1q1_ratio)[2], lty = 2)
      abline(h = get_median_CI(all_parameter_df$p0)[1], lty = 2)
      abline(h = get_median_CI(all_parameter_df$p0)[2], lty = 2)
      dev.off()
    }


    if(n_parameters >= 4){
      png(paste0(plot_folder_name,'/',fn_name,'/010.png'), width = 2000, height = 2000, res = 300)
      plot(parameter_df$r2q2_ratio, parameter_df$Psi_RTCDF,pch = 16, col = rgb(0,0,0,alpha = 0.1),
           main = paste0(title,'\n',fn_name,'\nr2/q2 vs Psi_RTCDF'),
           xlab = 'r2/q2',
           ylab = 'Psi_RTCDF')
      abline(v = get_CI(parameter_df$r2q2_ratio, 0.05)[1])
      abline(v = get_CI(parameter_df$r2q2_ratio, 0.05)[2])
      abline(h = get_CI(parameter_df$Psi_RTCDF, 0.05)[1])
      abline(h = get_CI(parameter_df$Psi_RTCDF, 0.05)[2])
      abline(v = get_median_CI(parameter_df$r2q2_ratio)[1], lty = 2)
      abline(v = get_median_CI(parameter_df$r2q2_ratio)[2], lty = 2)
      abline(h = get_median_CI(parameter_df$Psi_RTCDF)[1], lty = 2)
      abline(h = get_median_CI(parameter_df$Psi_RTCDF)[2], lty = 2)
      dev.off()


      png(paste0(plot_folder_name,'/',fn_name,'/011.png'), width = 2000, height = 2000, res = 300)
      plot(parameter_df$r2q2_ratio, parameter_df$p0,pch = 16, col = rgb(0,0,0,alpha = 0.1),
           main = paste0(title,'\n',fn_name,'\nr2q2 vs p0'),
           xlab = 'r2q2',
           ylab = 'p0')
      abline(v = get_CI(parameter_df$r2q2_ratio, 0.05)[1])
      abline(v = get_CI(parameter_df$r2q2_ratio, 0.05)[2])
      abline(h = get_CI(parameter_df$p0, 0.05)[1])
      abline(h = get_CI(parameter_df$p0, 0.05)[2])
      abline(v = get_median_CI(parameter_df$r2q2_ratio)[1], lty = 2)
      abline(v = get_median_CI(parameter_df$r2q2_ratio)[2], lty = 2)
      abline(h = get_median_CI(parameter_df$p0)[1], lty = 2)
      abline(h = get_median_CI(parameter_df$p0)[2], lty = 2)
      dev.off()

      png(paste0(plot_folder_name,'/',fn_name,'/012.png'), width = 2000, height = 2000, res = 300)
      plot(parameter_df$r1q1_ratio, parameter_df$r2q2_ratio,pch = 16, col = rgb(0,0,0,alpha = 0.1),
           main = paste0(title,'\n',fn_name,'\nr1q1 vs r2q2'),
           xlab = 'r1q1',
           ylab = 'r2q2')
      abline(v = get_CI(parameter_df$r1q1_ratio, 0.05)[1])
      abline(v = get_CI(parameter_df$r1q1_ratio, 0.05)[2])
      abline(h = get_CI(parameter_df$r2q2_ratio, 0.05)[1])
      abline(h = get_CI(parameter_df$r2q2_ratio, 0.05)[2])
      abline(v = get_median_CI(parameter_df$r1q1_ratio)[1], lty = 2)
      abline(v = get_median_CI(parameter_df$r1q1_ratio)[2], lty = 2)
      abline(h = get_median_CI(parameter_df$r2q2_ratio)[1], lty = 2)
      abline(h = get_median_CI(parameter_df$r2q2_ratio)[2], lty = 2)
      dev.off()
    }

  }
}

#' Write Parameter Table Helper Function
#'
#' This function generates table of optimized parameters
#' @param parameter_df Data frame of optimized parameters and other model function values (p0, Psi, etc)
#' @param folder_name Character vector indicating folder or directory name to be used when outputting table
#' @param model_fn_name Character vector used to indicate name of model function used for optimization
#' @param RGHD_m Int indicating m value of 2m-RGHD function if applicable
#' @export
write_parameter_table <- function(parameter_df, folder_name, model_fn_name, RGHD_m = 0){

  if(!(model_fn_name == 'RGHD')){
    fn_name <- str_replace(model_fn_name,'_',' ')

  }else{
    fn_name <- paste0('2m-RGHD (m=', RGHD_m, ')')
  }

  dir.create(paste0(folder_name,'/',fn_name))
  write.table(parameter_df, paste0(folder_name,'/',fn_name,'/',fn_name,' optimized_parameters.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}

#' Write Summary Table Helper Function
#'
#' This function generates summary statistics table of optimized parameters
#' @param parameter_df Data frame of optimized parameters and other model function values (p0, Psi, etc)
#' @param folder_name Character vector indicating folder or directory name to be used when outputting table
#' @param model_fn_name Character vector used to indicate name of model function used for optimization
#' @param RGHD_m Int indicating m value of 2m-RGHD function if applicable
#' @export
write_summary_table <- function(parameter_df, folder_name, model_fn_name, RGHD_m = 0){

  if(!(model_fn_name == 'RGHD')){
    fn_name <- str_replace(model_fn_name,'_',' ')

  }else{
    fn_name <- paste0('2m-RGHD (m=', RGHD_m, ')')
  }

  dir.create(paste0(folder_name,'/',fn_name))

  summary_df <- colnames(parameter_df) %>% as.data.frame()
  colnames(summary_df) <- 'Variable'
  summary_df$n <- nrow(parameter_df)
  summary_df$mean <- colMeans(parameter_df)
  summary_df$median <- colMedians(parameter_df %>% as.matrix())
  summary_df$std <- colSds(parameter_df %>% as.matrix())
  summary_df$mean_error_alpha_0.05 <- qt(1-0.05, summary_df$n-1) * summary_df$std / sqrt(summary_df$n)
  summary_df$mean_CI_alpha_0.05 <- paste0((summary_df$mean - summary_df$mean_error_alpha_0.05) %>% round(5), ', ',(summary_df$mean + summary_df$mean_error_alpha_0.05) %>% round(5))
  summary_df$median_CI_lower_rank <- ((summary_df$n/2) - (1.96 * sqrt(summary_df$n) / 2)) %>% round(0)
  summary_df$median_CI_upper_rank <- (1+ (summary_df$n/2) + (1.96 * sqrt(summary_df$n) / 2)) %>% round(0)
  summary_df$median_CI_approx_0.95 <- ''
  for(i in 1:length(parameter_df)){
    tmp <- sort(parameter_df[i] %>% unlist())
    summary_df$median_CI_approx_0.95[i] <- paste0(tmp[ summary_df$median_CI_lower_rank[i]]%>%round(5),', ',tmp[ summary_df$median_CI_upper_rank[i]]%>%round(5))
  }

  write.table(summary_df, paste0(folder_name,'/',fn_name,'/',fn_name,' parameter_summary.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}

#' Write Input Table Helper Function
#'
#' This function generates table of input data
#' @param folder_name Character vector indicating folder or directory name to be used when outputting table
#' @param data Vector of observed values
#' @export
write_input_table <- function(folder_name, data){

  input_df <- 1:length(data) %>% as.data.frame()
  colnames(input_df) <- 'x'
  input_df$data <- data

  write.table(input_df, paste0(folder_name,'/input_data.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}

#' SkeweDF Auto Helper Function
#'
#' This function will automatically optimize parameters for an empirical dataset given a model function and generate plots and tables
#' @param title Character vector indicating title of the empirical dataset, this will be present on every plot, this also determines the name of the folder where plots will be
#' @param data Vector of observed values
#' @param xlab Character vector indicating x axis label of plots, indicates what the random variable is
#' @param param_bounds A list of sequences which indicate space where parameters should be generated and fit
#' @param model_fn_name Character vector used to indicate name of model function used for optimization
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param n_cores Integer used to indicate number of cores to be used for this function if a socket cluster object is not defined.
#' @export
skeweDF_auto <- function(title = 'Dataset', data, xlab = 'Random Variable', param_bounds, model_fn_name, left_trunc = 1, right_trunc = left_trunc + length(data) - 1, n_cores = 1){
  dir.create('output')

  folder_name <- paste0('output/',title) # folder name processing
  # create all folders for all types of methods
  dir.create(folder_name)
  write_input_table(folder_name, data)

  if(model_fn_name == 'Kolmogorov Waring'){
    parameters <- local_fit_function(param_bounds = param_bounds, data = data[left_trunc:right_trunc], model_fn_name = 'Kolmogorov_Waring', n_cores, left_trunc = left_trunc, right_trunc = right_trunc)
    plot_model(title = title, model_fn_name = 'Kolmogorov_Waring',data = data, parameter_df = parameters[parameters$theta != 1,], n_parameters = 3, plot_folder_name = folder_name, xlab = xlab, left_trunc = left_trunc)
    write_parameter_table(parameter_df = parameters[parameters$theta != 1,],folder_name = folder_name,model_fn_name = 'Kolmogorov_Waring')
    write_summary_table(parameter_df = parameters[parameters$theta != 1,],folder_name = folder_name,model_fn_name = 'Kolmogorov_Waring')
  }
  else if(model_fn_name == 'RGHD'){
    m <- length(param_bounds) / 2
    RGHD_parameters <- local_fit_RGHD_ratio(param_bounds, data, n_cores, left_trunc = left_trunc, right_trunc = right_trunc)
    plot_model(title = title, model_fn_name = 'RGHD',data = data, parameter_df = RGHD_parameters, n_parameters = m*2, plot_folder_name = folder_name, xlab = xlab, left_trunc = left_trunc)
    write_parameter_table(parameter_df = RGHD_parameters,folder_name = folder_name,model_fn_name = 'RGHD',RGHD_m = m)
    write_summary_table(parameter_df = RGHD_parameters,folder_name = folder_name,model_fn_name = 'RGHD',RGHD_m = m)
  }
  else if(model_fn_name == 'Yule'){
    parameters <- local_fit_function(param_bounds = param_bounds, data = data[left_trunc:right_trunc], model_fn_name = 'Yule', n_cores, left_trunc = left_trunc, right_trunc = right_trunc)
    plot_model(title = title, model_fn_name = 'Yule',data = data, parameter_df = parameters, n_parameters = 1, plot_folder_name = folder_name, xlab = xlab, left_trunc = left_trunc)
    write_parameter_table(parameter_df = parameters,folder_name = folder_name,model_fn_name = 'Yule')
    write_summary_table(parameter_df = parameters,folder_name = folder_name,model_fn_name = 'Yule')
  }
  else if(model_fn_name == 'Generalized Yule'){
    parameters <- local_fit_function(param_bounds = param_bounds, data = data[left_trunc:right_trunc], model_fn_name = 'Generalized Yule', n_cores, left_trunc = left_trunc, right_trunc = right_trunc)
    plot_model(title = title, model_fn_name = 'Generalized Yule',data = data, parameter_df = parameters, n_parameters = 2, plot_folder_name = folder_name, xlab = xlab, left_trunc = left_trunc)
    write_parameter_table(parameter_df = parameters,folder_name = folder_name,model_fn_name = 'Generalized Yule')
    write_summary_table(parameter_df = parameters,folder_name = folder_name,model_fn_name = 'Generalized Yule')
  }
}

#' Psi Criterion Generalized KW
#'
#' This function calculates the Psi criterion goodness of fit metric given a set of parameters for the Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param pmf_weight Numeric of weight given to probability mass function for generation of Psi Criterion. For example, if pmf_weight <- 0.5, 50 percent of the Psi Criterion value will be attributed to the probability mass function while the other 50 percent will be attributed to the right-tail cumulative distribution function.
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @param ratio Boolean used to determine if parameters provide are a/b ratios rather than absolute values.
#' @export
psi_criterion_generalized_kw <- function(x, d, data, left_trunc, right_trunc, pmf_weight, weighted_rt, ratio = FALSE){
  if(ratio){
    x[1:d] <- x[1:d] * x[(d+1):(d+d)]
    print(x)
  }
  
  pmf <- psi_criterion_generalized_kw_pmf(x, d, data, left_trunc, right_trunc)
  cdf <- psi_criterion_generalized_kw_cdf(x, d, data, left_trunc, right_trunc, weighted_rt)

  return((pmf_weight * pmf) + ((1-pmf_weight)*cdf))
}


#' Psi Criterion Generalized KW PMF
#'
#' This function calculates the Psi criterion goodness of cit metric given a set of parameters for the probability mass function of the Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @export
psi_criterion_generalized_kw_pmf <- function(x, d, data, left_trunc, right_trunc){
  a <- x[1:d]
  b <- x[(d+1):(length(x)-1)]
  theta <- x[length(x)]

  model <- Kolmogorov_Waring(right_trunc, a, b, theta)
  model <- model[-1]
  model <- model[left_trunc:right_trunc]
  model <- model / sum(model)
  model <- model * sum(data)


  return(psi_criterion(data, model, length(x)) * -1)
}

#' Psi Criterion Generalized KW CDF
#'
#' This function calculates the Psi criterion goodness of cit metric given a set of parameters for the cumulative distribution function of the Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @export
psi_criterion_generalized_kw_cdf <- function(x, d, data, left_trunc, right_trunc, weighted_rt){
  a <- x[1:d]
  b <- x[(d+1):(length(x)-1)]
  theta <- x[length(x)]

  model <- Kolmogorov_Waring(right_trunc, a, b, theta)
  model <- model[-1]
  model <- model[left_trunc:right_trunc]
  model <- model / sum(model)
  model <- model * sum(data)

  if(weighted_rt){
    data <- weighted_right_tail_cdf(data)
    model <- weighted_right_tail_cdf(model)
  }
  else{
    data <- right_tail_cdf(data)
    model <- right_tail_cdf(model)
  }


  return(psi_criterion(data, model, length(x)) * -1)
}

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

#' Data Table to Sample
#' 
#' Helper function which converts a frequency dataframe to a sample vector
#' @param dt DataFrame. Frequency dataframe 
#'
#' @export
data_table_to_sample <- function(dt){
  out <- c()
  for(i in 1:nrow(dt)){
    out <- c(out,rep(dt[i,1], dt[i,2]))
  }
  return(out)
}

#' Sample to Data Table
#' 
#' Helper function which converts a sample vector to a frequency dataframe.
#' @param sample Vector
#'
#' @export
sample_to_data_table <- function(sample){
  out <- data.frame(table(sample))
  colnames(out) <- c('x','Freq')
  out[,1] <- out[,1] %>% as.character() %>% as.numeric()
  tmp <- 1:max(out[,1]) %>% as.data.frame()
  colnames(tmp) <- 'x'
  out <- left_join(tmp,out)
  out$Freq[is.na(out$Freq)] <- 0
  return(out)
}

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

#' Psi Criterion Mixed KW
#'
#' This function calculates the Psi criterion goodness of fit metric given a set of parameters for a mixed KW function
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @export
psi_criterion_mixed_kw <- function(x, d, data, left_trunc, right_trunc, pmf_weight, weighted_rt){
  pmf <- psi_criterion_mixed_kw_pmf(x, d, data, left_trunc, right_trunc)
  cdf <- psi_criterion_mixed_kw_cdf(x, d, data, left_trunc, right_trunc, weighted_rt)
  
  return((pmf_weight * pmf) + ((1-pmf_weight)*cdf))
}

#' Psi Criterion Mixed KW PMF
#'
#' This function calculates the Psi criterion goodness of fit metric given a set of parameters for the probability mass distribution function of the a mixed Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @export
psi_criterion_mixed_kw_pmf <- function(x, d, data, left_trunc, right_trunc){
  
  a1 <- x[1:d[1]]
  b1 <- x[(d[1]+1):(2*d[1])]
  theta1 <- x[(2*d[1]+1)]
  
  a2 <- x[(2*d[1]+2):(2*d[1]+2+d[2]-1)]
  b2 <- x[(2*d[1]+2+d[2]):(2*d[1]+2+2*d[2]-1)]
  theta2 <- x[2*d[1]+2+2*d[2]]
  
  model1_weight <- unlist(x[2*d[1]+2+2*d[2]+1])
  
  model1 <- Kolmogorov_Waring(right_trunc, unlist(a1), unlist(b1), unlist(theta1))
  model2 <- Kolmogorov_Waring(right_trunc, unlist(a2), unlist(b2), unlist(theta2))
  model <- additive_mixed_model(model1, model2, model1_weight)
  
  model <- model[-1]
  model <- model[left_trunc:right_trunc]
  model <- model / sum(model)
  model <- model * sum(data)
  
  
  return(psi_criterion(data, model, length(x)) * -1)
}

#' Psi Criterion Mixed KW CDF
#'
#' This function calculates the Psi criterion goodness of fit metric given a set of parameters for the cumulative distribution function of the a mixed Generalized Kolmogorov Waring function.
#' @param x Vector of parameters
#' @param d Int used to indicate the number of 'a' terms within the 'param_bounds' variable. The remaining values will be considered 'b' terms.
#' @param data Vector of observed values
#' @param left_trunc Int used to determine starting index of model to use for optimization
#' @param right_trunc Int used to determine ending index of model to use for optimization
#' @param weighted_rt Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.
#' @export
psi_criterion_mixed_kw_cdf <- function(x, d, data, left_trunc, right_trunc, weighted_rt){
  a1 <- x[1:d[1]]
  b1 <- x[(d[1]+1):(2*d[1])]
  theta1 <- x[(2*d[1]+1)]
  
  a2 <- x[(2*d[1]+2):(2*d[1]+2+d[2]-1)]
  b2 <- x[(2*d[1]+2+d[2]):(2*d[1]+2+2*d[2]-1)]
  theta2 <- x[2*d[1]+2+2*d[2]]
  
  model1_weight <- unlist(x[2*d[1]+2+2*d[2]+1])
  
  model1 <- Kolmogorov_Waring(right_trunc, unlist(a1), unlist(b1), unlist(theta1))
  model2 <- Kolmogorov_Waring(right_trunc, unlist(a2), unlist(b2), unlist(theta2))
  model <- additive_mixed_model(model1, model2, model1_weight)
  model <- model[-1]
  model <- model[left_trunc:right_trunc]
  model <- model / sum(model)
  model <- model * sum(data)
  
  if(weighted_rt){
    data <- weighted_right_tail_cdf(data)
    model <- weighted_right_tail_cdf(model)
  }
  else{
    data <- right_tail_cdf(data)
    model <- right_tail_cdf(model)
  }
  
  
  return(psi_criterion(data, model, length(x)) * -1)
}

#' Fit Mixed Kolmogorov-Waring
#'
#'#' This function will take in a set of parameters and calculate goodness of fit based on the Psi criterion metric
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
#' @export
fit_mixed_kw <- function(param_bounds,
                         d,
                         data,
                         weighted_rt = FALSE,
                         pmf_weight = 0.0,
                         par_chunk = 100,
                         par_chunk_size = 10,
                         n_cores = 12,
                         clust,
                         left_trunc = 1,
                         right_trunc = left_trunc+length(data)-1){
  
  data <- data / sum(data)
  
  standalone <- FALSE
  
  
  if(missing(clust)){
    standalone <- TRUE
    clust <- makeCluster(n_cores)
    clusterExport(clust, varlist = c('data','param_bounds','weighted_rt', 'par_chunk_size','d','left_trunc','right_trunc', 'pmf_weight'), envir = environment())
    clusterCall(clust, function() library(SkeweDF))
    clusterCall(clust, function() library(tidyverse))
  }
  
  clusterExport(clust, varlist = c('data','param_bounds','weighted_rt', 'par_chunk_size','d','left_trunc','right_trunc', 'pmf_weight'), envir = environment())
  
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
    
    fn_parameters <- multistart(parmat = par_mat,fn = psi_criterion_mixed_kw,method = 'L-BFGS-B',
                                lower = param_lower, upper = param_upper,
                                data = data, d = d, pmf_weight = pmf_weight, weighted_rt = weighted_rt,
                                left_trunc = left_trunc, right_trunc = right_trunc);
    
    fn_parameters$fevals <- NULL
    fn_parameters$gevals <- NULL
    fn_parameters$convergence <- NULL
    
    for(i in 1:d[1]){
      colnames(fn_parameters)[i] <- paste0('a1_',i)
    }
    for(i in (d[1]+1):(2*d[1])){
      colnames(fn_parameters)[i] <- paste0('b1_',i-d)
    }
    colnames(fn_parameters)[2*d[1]+1] <- paste0('theta1')
    
    for(i in (2*d[1]+2):(2*d[1]+2+d[2]-1)){
      colnames(fn_parameters)[i] <- paste0('a1_',(i-(2*d[1]+d[2]-1)))
    }
    for(i in (2*d[1]+2+d[2]):(2*d[1]+2+2*d[2]-1)){
      colnames(fn_parameters)[i] <- paste0('b1_',(i-(2*d[1]+2*d[2]-1)))
    }
    colnames(fn_parameters)[2*d[1]+2+2*d[2]] <- paste0('theta2')
    #fn_parameters[1:m] <- t(apply(fn_parameters[1:m], 1, FUN=function(x) sort(x)))
    #fn_parameters[(m+1):(length(fn_parameters)-2)] <- t(apply(fn_parameters[(m+1):(length(fn_parameters)-2)], 1, FUN=function(x) sort(x)))
    
    colnames(fn_parameters)[length(fn_parameters)-1] <- 'model1_weight'
    colnames(fn_parameters)[length(fn_parameters)] <- 'Psi'
    fn_parameters$Psi <- fn_parameters$Psi * -1
    fn_parameters$pmf_weight <- pmf_weight
    
    return(fn_parameters)
  }) %>% bind_rows()
  
  print('Complete')
  
  if(standalone){
    stopCluster(clust)
  }
  
  return(parameters)
  
}

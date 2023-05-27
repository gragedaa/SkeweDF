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

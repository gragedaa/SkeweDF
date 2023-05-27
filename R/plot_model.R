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

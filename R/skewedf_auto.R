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

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

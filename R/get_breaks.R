get_breaks <- function(years){
  
  # Define bin width
  bin_width <- 5
  
  # Get the data range
  data_range <- range(years)
  
  # Calculate "nice" lower and upper bounds
  lower <- floor(data_range[1] / bin_width) * bin_width
  upper <- ceiling(data_range[2] / bin_width) * bin_width
  
  # Create breaks every 5 years
  breaks <- seq(lower, upper, by = bin_width)
  
  return(breaks)
}
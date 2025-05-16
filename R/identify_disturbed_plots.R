identify_disturbed_plots <- function(df){
  
  # identify disturbed plots based on whether plot of at least one difference 
  # between inventories looks like disturbance
  plots_disturbed <- df |> 
    group_by(plotID) |> 
    nest() |> 
    mutate(data = purrr::map(data, ~identify_disturbed_byinventory(.))) |> 
    mutate(ndisturbed = purrr::map_int(data, ~get_sum_disturbed(.))) |> 
    unnest(data)
}

get_sum_disturbed <- function(df){
  sum(df$disturbed, na.rm = TRUE)
}

# counts intervals between inventories that look like disturbance was at play
identify_disturbed_byinventory <- function(df){
  df |> 
    mutate(
      dlogQMD = logQMD - lag(logQMD),
      dlogDensity = logDensity - lag(logDensity)
    ) |> 
    mutate(
      disturbed = ifelse(dlogQMD < 0 & dlogDensity < 0, TRUE, FALSE)
      )
    
}

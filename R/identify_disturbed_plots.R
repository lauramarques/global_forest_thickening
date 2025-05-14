identify_disturbed_plots <- function(df){
  
  # identify disturbed plots based on whether plot of at least one difference 
  # between inventories looks like disturbance
  plots_disturbed <- df |> 
    group_by(plotID) |> 
    nest() |> 
    mutate(ndisturb = purrr::map_int(data, ~get_ndisturb(.))) |> 
    filter(ndisturb > 0) |> 
    pull(plotID) |> 
    unique()
  
  df <- df |> 
    mutate(disturbed = ifelse(plotID %in% plots_disturbed, TRUE, FALSE))
}

# counts intervals between inventories that look like disturbance was at play
get_ndisturb <- function(df){
  ndisturb <- df |> 
    dplyr::select(year, logQMD, logDensity) |> 
    mutate(
      dlogQMD = logQMD - lag(logQMD),
      dlogDensity = logDensity - lag(logDensity)
    ) |> 
    filter(dlogQMD < 0 & dlogDensity < 0) |> 
    nrow()
}

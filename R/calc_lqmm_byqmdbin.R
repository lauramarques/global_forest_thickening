calc_lqmm_byqmdbin <- function(df){
  
  # make sure at least 300 points are on average in each bin
  nbins <- max(round(length(df$logQMD[!is.na(df$logQMD)]) / 300), 10)
  bin_edges <- pretty(df$logQMD, n = nbins)
  bin_labels <- bin_edges[1:length(bin_edges)-1] + (bin_edges[2] - bin_edges[1])/2
  
  df |> 
    mutate(bin_lqmm = cut(
      logQMD, 
      breaks = bin_edges, 
      length.out = nbins + 1, 
      include.lowest = TRUE,
      labels = bin_labels
      )
    ) |> 
    group_by(bin_lqmm) |> 
    nest() |> 
    mutate(nvals = purrr::map_int(data, ~nrow(.))) |> 
    filter(nvals >= 30) |> 
    mutate(mod_lqmm = purrr::map(data, ~lqmm(
      logDensity ~ year,
      random = ~1,
      group = plotID,
      tau = c(0.90),
      type = "normal",
      data = .
    ))) |> 
    mutate(sum_lqmm = purrr::map(
      mod_lqmm,
      ~try(summary(.)))) |> 
    mutate(failed = purrr::map_lgl(
      sum_lqmm,
      ~determine_failed(.))) |> 
    filter(!failed) |> 
    mutate(pval = purrr::map_dbl(
      sum_lqmm,
      ~get_pval(.)
    )) |> 
    mutate(coef_year = purrr::map_dbl(
      sum_lqmm,
      ~get_coef_year(.)
    )) |> 
    mutate(coef_year_upper = purrr::map_dbl(
      sum_lqmm,
      ~get_coef_year_upper(.)
    )) |> 
    mutate(coef_year_lower = purrr::map_dbl(
      sum_lqmm,
      ~get_coef_year_lower(.)
    )) |> 
    mutate(upwardshift = ifelse(coef_year > 0 & pval < 0.01, TRUE, FALSE))
  
}

get_pval <- function(sum){
  sum$tTable["year", "Pr(>|t|)"]
}
get_coef_year <- function(sum){
  sum$tTable["year", "Value"]
}
get_coef_year_upper <- function(sum){
  sum$tTable["year", "upper bound"]
}
get_coef_year_lower <- function(sum){
  sum$tTable["year", "lower bound"]
}
determine_failed <- function(out){
  class(out) == "try-error"
}


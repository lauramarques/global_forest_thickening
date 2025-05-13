plot_lqmm_bybiome <- function(data, mod, name){
  
  df_pred <- purrr::map_dfr(
    as.list(c(1985, 2000, 2015)),
    ~ tibble(
      logQMD = seq(from = min(data$logQMD), to = max(data$logQMD), length.out = 100),
      year = .x,
      plotID = unique(data$plotID)[1]  # may take any value here, as long as predict(..., level = 0)
    )
  ) %>%
    mutate(
      pred_70 =  predict(
        mod,
        level = 0,
        newdata = .
      )[,"0.7"],
      pred_90 =  predict(
        mod,
        level = 0,
        newdata = .
      )[,"0.9"]
    ) |> 
    mutate(year = as.factor(year))
  
  ggplot() + 
    geom_point(
      data = data, 
      aes(x = logQMD, y = logDensity), 
      alpha = 0.7, 
      col = "darkgrey", 
      size = 1, 
      shape = 16, 
      inherit.aes = FALSE
    ) + 
    geom_line(
      data = df_pred,
      aes(logQMD, pred_70, color = year, group = year)
    ) +
    labs(
      x = expression(ln(QMD)),
      y = expression(ln(italic(N))), 
      title = name,
      subtitle = bquote(
        italic(n) == .(as.character(nobs(mod)))
      ),
      color  = "Year"
    ) +
    scale_color_manual("Year",
                       breaks = c(1985, 2000, 2015), 
                       values = c(viridis(3)[2], viridis(3)[3], viridis(3)[1])) +
    scale_fill_manual("Year",
                      breaks = c(1985, 2000, 2015), 
                      values = c(viridis(3)[2], viridis(3)[3], viridis(3)[1])) +
    theme_classic()
  
}

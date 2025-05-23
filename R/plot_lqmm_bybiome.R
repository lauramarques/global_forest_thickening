plot_lqmm_bybiome <- function(data, mod, name, plot_legend = FALSE){
  
  # predict STL as 70% and 90% quantiles from regression fit for three chosen years
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
  
  # extract fit info
  out <- summary(mod)
  df_caption <- tibble(
    variable = c("year"),
    estimate = out$tTable$`0.9`[c("year"), "Value"],
    pvalue = out$tTable$`0.9`[c("year"), "Pr(>|t|)"]
  ) |> 
    mutate(
      estimate_lab = round(estimate, 3),
      pval_lab = ifelse(pvalue < 0.001, "< 0.001 ***", paste0(round(pvalue, 3)))
    )
  
  # construct subtitle
  subtitle <- bquote(
    italic(n) == .(as.character(nobs(mod))) ~ ~ ~
      italic(p)(yr) == .(df_caption$pval_lab[1])
  )
  
  ggplot() + 
    geom_point(
      data = data, 
      aes(x = logQMD, y = logDensity), 
      col = "grey30", 
      alpha = 0.1,
      size = 1, 
      shape = 16, 
      inherit.aes = FALSE
    ) + 
    geom_line(
      data = df_pred,
      aes(logQMD, pred_90, color = year, group = year)
    ) +
    labs(
      x = expression(ln(QMD)),
      y = expression(ln(italic(N))), 
      title = name,
      subtitle = subtitle,
      color  = "Year"
    ) +
    scale_color_manual("Year",
                       breaks = c(1985, 2000, 2015), 
                       values = c(viridis(3)[2], viridis(3)[3], viridis(3)[1])) +
    scale_fill_manual("Year",
                      breaks = c(1985, 2000, 2015), 
                      values = c(viridis(3)[2], viridis(3)[3], viridis(3)[1])) +
    theme_classic() +  
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12), 
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10), 
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 11),
      legend.key = element_rect(fill = NA, color = NA),
      legend.position = ifelse(plot_legend, "bottom", "none"),
      plot.caption = element_text(vjust = -1),
      plot.title.position = "plot"
    )
    # scale_x_continuous(limits = c(2,4.5), breaks = seq(2,4,1))
    # scale_y_continuous(limits = c(2.9,9.3), breaks = seq(4,8,2))
  
}

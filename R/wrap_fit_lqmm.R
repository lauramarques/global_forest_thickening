wrap_fit_lqmm <- function(split){
  
  # Fit model — suppress warnings (e.g., convergence issues)
  tryCatch({
    mod <- lqmm(
      logDensity ~ logQMD + year,
      random = ~1,
      group = plotID,
      tau = 0.90,
      data = analysis(split),
      type = "normal",
      control = list(
        LP_max_iter = 500,
        LP_tol_ll = 1e-4
      )
    )
    
    tibble(term = names(coef(mod)), estimate = coef(mod))
  }, error = function(e) NULL)
}
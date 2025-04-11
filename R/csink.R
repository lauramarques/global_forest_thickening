csink <- function(data_all, fit, a_mean, a_sd){
  # Bootstrap with loop ...
  
  # 1. Sample QMD from its distribution in the data => QMDj
  # QMDj <- sample(data_all$QMD, 1)
  # logQMDj <- log(QMDj)
  # 1*. Sample the row and het QMD and the random variables
  random_row_index <- sample(1:nrow(data_all), 1)
  random_values <- data_all[random_row_index, c("QMD", "species", "dataset", "plotID")]
  QMDj <- random_values$QMD
  
  # 2. Estimate mean N given QMDj using the LMM  (lnN ~ lnQMD + year + 1|...) and two years (for example year = 2000, 2001) => N0, N1.
  # N0_mean <- as.data.frame(ggpredict(fit1, terms = c("logQMD[logQMDj]","year[2000]"), full.data = TRUE))$predicted
  # N1_mean <- as.data.frame(ggpredict(fit1, terms = c("logQMD[logQMDj]","year[2001]"), full.data = TRUE))$predicted
  # Single new data point t0
  point_t0 <- random_values |> 
    mutate(logQMD = log(QMDj),
           year = 2000)
  
  # Predict for the new point, including random effects
  logN0_mean <- predict(fit, newdata = point_t0, re.form = NULL, allow.new.levels = TRUE) # re.form = NULL includes random effects, while re.form = NA removes random effects, like ggpredict(..., full.data = FALSE).
  N0_mean <- exp(logN0_mean)
  
  # Single new data point t1
  point_t1 <- point_t0 |> mutate(year = 2001)
  
  # Predict for the new point, including random effects
  logN1_mean <- predict(fit, newdata = point_t1, re.form = NULL, allow.new.levels = TRUE)
  N1_mean <- exp(logN1_mean)
  
  # 3. Estimate the biomass change, given the QMDj, N0 and N1 as: dB = ak * QMDj^2 * (N1 - N0), 
  # where ak is a sampled value from the fitted a (a_mean), considering its standard error (a_sd) and a normal distribution.
  ak <- rnorm(1, a_mean, a_sd)
  dB = ak * QMDj^2 * (N1_mean - N0_mean)
  
  out <- tibble(QMDj = QMDj, 
                N0 = N0_mean, 
                N1 = N1_mean, 
                ak = ak, 
                dB_Mg_ha = dB * 10^-3) 
  
  return(out)
}

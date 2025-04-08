csink_global <- function(
    global_drivers_row, data_all,
    coef_ai_mean, coef_ai_sd,
    coef_ndep_mean, coef_ndep_sd,
    coef_orgc_mean, coef_orgc_sd,
    coef_pbr_mean, coef_pbr_sd,
    coef_aiyear_mean, coef_aiyear_sd,
    coef_ndepyear_mean, coef_ndepyear_sd,
    coef_orgcyear_mean, coef_orgcyear_sd,
    coef_pbryear_mean, coef_pbryear_sd){
  
  # 1*. Sample the row and het QMD and the random variables
  random_row_index <- sample(1:nrow(data_all), 1)
  random_values <- data_all[random_row_index, c("QMD", "species", "dataset", "plotID")]
  QMDj <- random_values$QMD
  
  # Sample from the model coefficients for MI, Ndep, ORGC, PBR
  b_ai <- rnorm(1, coef_ai_mean, coef_ai_sd)
  b_ndep <- rnorm(1, coef_ndep_mean, coef_ndep_sd)
  b_orgc <- rnorm(1, coef_orgc_mean, coef_orgc_sd)
  b_pbr <- rnorm(1, coef_pbr_mean, coef_pbr_sd)
  b_aiyear <- rnorm(1, coef_aiyear_mean, coef_aiyear_sd)
  b_ndepyear <- rnorm(1, coef_ndepyear_mean, coef_ndepyear_sd)
  b_orgcyear <- rnorm(1, coef_orgcyear_mean, coef_orgcyear_sd)
  b_pbryear <- rnorm(1, coef_pbryear_mean, coef_pbryear_sd)
  
  # 2. Estimate mean N given QMDj using the LMM  (lnN ~ lnQMD + year + 1|...) and two years (for example year = 2000, 2001) => N0, N1.
  # Single new data points for t0 and t1
  
  # Important: Need to manually scale new predictors using the same means and standard deviations as in data_all
  # Compute means and standard deviations from the original data
  means <- colMeans(data_all[, c("logQMD", "year", "ai", "ndep", "ORGC", "PBR")], na.rm = TRUE)
  sds <- apply(data_all[, c("logQMD", "year", "ai", "ndep", "ORGC", "PBR")], 2, sd, na.rm = TRUE)
  
  # Create newdata with a single row
  point_t0 <- data.frame(
    logQMD = (log(QMDj) - means["logQMD"]) / sds["logQMD"],  
    year   = (2000 - means["year"]) / sds["year"],  
    ai     = (global_drivers_row$ai - means["ai"]) / sds["ai"],
    ndep   = (global_drivers_row$ndep - means["ndep"]) / sds["ndep"],
    ORGC   = (global_drivers_row$ORGC - means["ORGC"]) / sds["ORGC"],
    PBR    = (global_drivers_row$PBR - means["PBR"]) / sds["PBR"]
  )
  
  point_t1 <- point_t0 |> 
    mutate(year = (2001 - means["year"]) / sds["year"])
  
  # Extract design matrix (without scale() )
  X <- model.matrix(~ logQMD + year * ai + year * ndep + year * ORGC + year * PBR, 
                    point_t0)
  #print(X)
  
  Y <- model.matrix(~ logQMD + year * ai + year * ndep + year * ORGC + year * PBR, 
                    point_t1)
  #print(Y)
  
  # Extract original fixed effects
  # original_fixed_effects <- fixef(fit1)
  
  # Use modified fixed effect coefficients
  new_coeffs <- c(coef_intercept,
                  coef_logQMD,
                  coef_year,
                  b_ai,
                  b_ndep,
                  b_orgc,
                  b_pbr,
                  b_aiyear,
                  b_ndepyear,
                  b_orgcyear,
                  b_pbryear)  
  
  # Ensure dimensions match
  #print(length(new_coeffs))  # Should match ncol(X)
  #print(dim(X))  # Should be (1, number_of_predictors + 1)
  
  # Compute predictions for point_t0 and new sampled coefficients
  # logN0_mean <- predict(fit1, newdata = point_t0, re.form = NULL, allow.new.levels = TRUE)
  logN0_mean <- X %*% new_coeffs
  N0_mean <- exp(logN0_mean)
  
  # Compute predictions for point_t1 and new sampled coefficients
  #logN1_mean <- predict(fit1, newdata = point_t1, re.form = NULL, allow.new.levels = TRUE)
  logN1_mean <- Y %*% new_coeffs
  N1_mean <- exp(logN1_mean)
  
  # 3. Estimate the biomass change, given the QMDj, N0 and N1 as: dB = ak * QMDj^2 * (N1 - N0), 
  # where ak is a sampled value from the fitted a (a_mean), considering its standard error (a_sd) and a normal distribution.
  ak <- rnorm(1, a_mean, a_sd)
  dB = ak * QMDj^2 * (N1_mean - N0_mean)
  
  out <- tibble(QMDj = QMDj, 
                b_ai = b_ai, 
                b_ndep = b_ndep,
                b_orgc = b_orgc,
                b_pbr = b_pbr,
                b_aiyear = b_aiyear,
                b_ndepyear = b_ndepyear,
                b_orgcyear = b_orgcyear,
                b_pbryear = b_pbryear,
                N0 = N0_mean, 
                N1 = N1_mean, 
                a = ak, 
                dB_Mg_ha = dB * 10^-3,
                lon = global_drivers_row$lon, 
                lat = global_drivers_row$lat, 
                area_ha = global_drivers_row$area_ha)
  
  return(out)
}
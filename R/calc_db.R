# Everything below is for one row and should be done rowwise -------
calc_db <- function(
    df_samples_onerow,
    data_forest_plots_selfthinning_means,
    data_forest_plots_selfthinning_sds,
    coef_samples_selfthinning,
    coef_samples_biomass
    ){
  
  # Create newdata with a single row and scale with mean and SD from data used 
  # for fitting the selfthinning model
  # xxx test
  tmp_t0 <- df_samples_onerow |> 
    dplyr::mutate(year = 2000) |> 
    dplyr::select(logQMD, year, ai, ndep, ORGC, PBR)
  
  tmp_t1 <- df_samples_onerow |> 
    dplyr::mutate(year = 2001) |> 
    dplyr::select(logQMD, year, ai, ndep, ORGC, PBR)
  
  # normalise
  tmp_t0_norm <- (tmp_t0 - data_forest_plots_selfthinning_means)/data_forest_plots_selfthinning_sds
  tmp_t1_norm <- (tmp_t1 - data_forest_plots_selfthinning_means)/data_forest_plots_selfthinning_sds
  
  # Extract design matrix XXX why without scale()?
  X <- model.matrix(
    ~ logQMD + year * ai + year * ndep + year * ORGC + year * PBR, 
    tmp_t0_norm
  )
  Y <- model.matrix(
    ~ logQMD + year * ai + year * ndep + year * ORGC + year * PBR, 
    tmp_t1_norm
  )
  
  new_coeffs_selfthinning <- df_samples_onerow |> 
    dplyr::select(all_of(names(as.data.frame(coef_samples_selfthinning)))) |> 
    as.matrix()
  
  # # Ensure dimensions match
  # print(length(new_coeffs_selfthinning))  # Should match ncol(X)
  # print(dim(X))  # Should be (1, number_of_predictors + 1)
  
  # Compute predictions for point_t0 and new sampled coefficients
  logN0_mean <- X %*% as.numeric(new_coeffs_selfthinning)
  N0_mean <- c(exp(logN0_mean))
  
  # Compute predictions for point_t1 and new sampled coefficients
  #logN1_mean <- predict(fit1, newdata = point_t1, re.form = NULL, allow.new.levels = TRUE)
  logN1_mean <- Y %*% as.numeric(new_coeffs_selfthinning)
  N1_mean <- c(exp(logN1_mean))
  
  # 3. Estimate the biomass change, given the QMDj, N0 and N1 as: dB = ak * QMDj^2 * (N1 - N0), 
  ak <- df_samples_onerow |> 
    dplyr::select(all_of(names(as.data.frame(coef_samples_biomass)))) |> 
    as.matrix() |> 
    c()
  
  qmd_j <- df_samples_onerow |> 
    dplyr::pull(logQMD) |> 
    exp()
  
  db <- ak * qmd_j^2 * (N1_mean - N0_mean)
  
  db <- ifelse(length(db) == 0, NA, db)
  
  return(db)
}

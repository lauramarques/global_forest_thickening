library(tidyverse)

data_fil_biome <- readRDS(here::here("data/data_fil_biomes.rds"))

data_fil_biome <- data_fil_biome |>
  mutate(NQMD2 = density * QMD^2)  # new variable for modelling biomass 

## Biome 1 ---------
# Tropical & Subtropical Moist Broadleaf Forests
df_sub <- data_fil_biome |> 
  filter(biomeID == 1)

fit_selfthinning = lmer(
  log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
  data = df_sub, 
  na.action = "na.exclude"
)

fit_biomass = lmer(
  biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
  data = df_sub, 
  na.action = "na.exclude"
)

# Single iteration - pack into a function
# Sample one row
set.seed(1982)
df_sub_sample <- df_sub |> 
  slice_sample(n = 1)
  
# Estimate mean N given QMDj using the LMM  (lnN ~ lnQMD + year + 1|...) and two years (for example year = 2000, 2001) => N0, N1.
# New data points for two years (1 year apart)
newdata_t0 <- data.frame(logQMD = df_sub_sample$logQMD, year = 2000, dataset = "dummy", plotID = "dummy", species = "dummy")
newdata_t1 <- newdata_t0 |> 
  mutate(year = 2001)

# Predict the fixed effects part only (exclude random effects). 
# (marginal) predictions based only on fixed effects.
log_density_t0 <- predict(fit_selfthinning, newdata = newdata_t0, re.form = NA)
log_density_t1 <- predict(fit_selfthinning, newdata = newdata_t1, re.form = NA)

# back-transform
density_t0 <- exp(log_density_t0)
density_t1 <- exp(log_density_t1)

# Estimate the biomass change, given the QMDj, N0 and N1 as: dB = ak * QMDj^2 * (N1 - N0), 
# where ak is a sampled value from the fitted a (a_mean), considering its standard error (a_sd) and a normal distribution.
newdata_biomass <- data.frame(
  density = c(density_t0, density_t1),
  QMD = rep(df_sub_sample$QMD, 2),
  dataset = rep("dummy", 2), 
  plotID = rep("dummy", 2)
  ) |> 
  mutate(NQMD2 = density * QMD^2) %>%   # new variable for modelling biomass 
  mutate(biomass_pred_mean = predict(fit_biomass, newdata = ., re.form = NA))  # (marginal) predictions based only on fixed effects.

# new version: considering residual error (prediction error) in biomass estimate, given N*QMD^2, 
# and assuming uncorrelated errors in the two biomass estimates, not considering random effect uncertainty 
dB <- newdata_biomass$biomass_pred_mean[2] - newdata_biomass$biomass_pred_mean[1]
dB_with_unc <- dB + rnorm(1, mean = 0, sd = sqrt(2) * sigma(fit_biomass))

# # original version: considering uncertainty in coefficient relating N*QMD^2 to biomass.
# a_mean <- summary(fit_biomass)$coefficient[1,1]
# a_sd <- summary(fit_biomass)$coefficient[1,2]
# ak <- rnorm(1, a_mean, a_sd)
# dB_with_unc <- ak * newdata_biomass$QMD^2 * (newdata_biomass$density[2] - newdata_biomass$density[1])


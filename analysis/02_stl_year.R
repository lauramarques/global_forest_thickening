# This script analyses the changes in the STLs over time (calendar year) by biome.

# Load packages ----------------------------------------------------------------
#library(renv)
library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(rFIA)
library(patchwork)
library(terra)
library(sf)
library(lme4)
library(lmerTest)
library(ggeffects)
library(effects)
library(sjPlot)
library(measurements)
library(sp)
library(lqmm)
library(ggforce)
library(MuMIn)
# library(ingestr)
library(DescTools)
library(here)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

# Load functions ---------------------------------------------------------------
source(here("R/functions.R"))
source(here("R/plot_stl_bybiome.R"))

# load data
data_fil_biomes <- readRDS(here("data/data_fil_biomes.rds"))
data_fil_biomes <- readRDS(here("~/GFDYglobe/data/inputs/data_fil_biomes.rds"))

plot_map_fil <-  plot_map(data_fil_biomes)
plot_map_fil

# get temporally fixed self-thinning line
plot_stl_fil <- plot_stl(data_fil_biomes)
plot_stl_fil

# Biome 1: Tropical & Subtropical Moist Broadleaf Forests  ---------------------

# XXX In my interpretation because I don't have the file above, and as far as I 
# can see, these are identical:
# waldo::compare(
#   readRDS(here::here("data/data_fil_biome1.rds")), 
#   data_fil_biomes |>
#     filter(biomeID == 1)
#   )

data_fil_biome1 <- data_fil_biomes |>
  filter(biomeID == 1)

## Linear mixed effects model --------------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_biome1 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome1, 
  na.action = "na.exclude"
  )
summary(mod_lmm_biome1)

### Plot self-thinning line ----------------------------------------------------
gg_stl_biome1 <- plot_stl_bybiome(
  data_fil_biome1, 
  mod_lmm_biome1, 
  name = "Tropical Moist Broadleaf Forests", 
  years = c(1985, 2000, 2015)
) + theme(legend.position = "right")

### Data over years ------------------------------------------------------------
hist_Year <- ggplot(data_fil_biome1, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic()
hist_Year

### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome1, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
  )

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
    ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome1)
# print(out)
# r.squaredGLMM(mod_lmm_biome1)
# 
# plot(allEffects(mod_lmm_biome1))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome1,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
#   ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome1$year))
# years


# Biome 4: Temperate Broadleaf & Mixed Forests  --------------------------------

data_fil_biome4 <- data_fil_biomes |>
  filter(biomeID == 4)

## Linear mixed effects model --------------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_biome4 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome4, 
  na.action = "na.exclude"
)
summary(mod_lmm_biome4)

### Plot self-thinning line ----------------------------------------------------
gg_stl_biome4 <- plot_stl_bybiome(
  data_fil_biome4, 
  mod_lmm_biome4, 
  name = "Temperate Broadleaf & Mixed Forests", 
  years = c(1985, 2000, 2015)
)

### Data over years ------------------------------------------------------------
hist_Year <- ggplot(data_fil_biome4, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic()

hist_Year

### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome4, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome4)
# print(out)
# r.squaredGLMM(mod_lmm_biome4)
# 
# plot(allEffects(mod_lmm_biome4))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome4,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome4$year))
# years

# Biome 5: Temperate Conifer Forests Forest  --------------------------------
data_fil_biome5 <- data_fil_biomes |>
  filter(biomeID == 5)

## Linear mixed effects model --------------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_biome5 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome5, 
  na.action = "na.exclude"
)
summary(mod_lmm_biome5)

### Plot self-thinning line ----------------------------------------------------
gg_stl_biome5 <- plot_stl_bybiome(
  data_fil_biome5, 
  mod_lmm_biome5, 
  name = "Temperate Conifer Forests Forest", 
  years = c(1985, 2000, 2015)
)

### Data over years ------------------------------------------------------------
hist_Year <- ggplot(data_fil_biome5, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic()

hist_Year

### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome5, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome5)
# print(out)
# r.squaredGLMM(mod_lmm_biome5)
# 
# plot(allEffects(mod_lmm_biome5))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome5,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome2$year))
# years

# Biome 6: Boreal Forests/Taiga Forest  --------------------------------
data_fil_biome6 <- data_fil_biomes |>
  filter(biomeID == 6)

## Linear mixed effects model --------------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_biome6 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome6, 
  na.action = "na.exclude"
)

### Plot self-thinning line ----------------------------------------------------
gg_stl_biome6 <- plot_stl_bybiome(
  data_fil_biome6, 
  mod_lmm_biome6, 
  name = "Boreal Forests/Taiga Forest", 
  years = c(1985, 2000, 2015)
)

### Data over years ------------------------------------------------------------
hist_Year <- ggplot(data_fil_biome6, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic()

hist_Year

### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome6, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome6)
# print(out)
# r.squaredGLMM(mod_lmm_biome6)
# 
# plot(allEffects(mod_lmm_biome6))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome6,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome2$year))
# years


# Biome 12: Mediterranean Forests ----------------------------------------------
data_fil_biome12 <- data_fil_biomes |>
  filter(biomeID == 12)

## Linear mixed effects model --------------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_biome12 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome12, 
  na.action = "na.exclude"
)

### Plot self-thinning line ----------------------------------------------------
gg_stl_biome12 <- plot_stl_bybiome(
  data_fil_biome12, 
  mod_lmm_biome12, 
  name = "Mediterranean Forests", 
  years = c(1985, 2000, 2015)
)

### Data over years ------------------------------------------------------------
hist_Year <- ggplot(data_fil_biome12, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic()

hist_Year

### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome12, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome12)
# print(out)
# r.squaredGLMM(mod_lmm_biome12)
# 
# plot(allEffects(mod_lmm_biome12))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome12,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome2$year))
# years
# 

# Publication Figure 1 ---------------------------------------------------------
fig1 <- 
  gg_stl_biome1 + 
  gg_stl_biome4 + 
  gg_stl_biome5 + 
  gg_stl_biome6 + 
  gg_stl_biome12 +
  guide_area() +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "a",tag_suffix = ")")

ggsave(
  filename = here::here("manuscript/figures/fig1.png"),
  plot = fig1,
  width = 11, 
  height = 7.5
)

# # Special case of BCI ----------------------------------------------------------
# 
# data_bci <- readRDS(here::here("data/inputs/data_bci.rds"))
# data_bci
# plot_stl(data_bci)
# plot_map(data_bci)
# 
# data_bci_1_3 <- data_bci |> 
#   filter(census <= 3)
# 
# data_bci_1_3 <- data_bci |> 
#   filter(census ==1 |
#          census ==2|
#          census ==6)
# 
# 
# data_bci_4_8 <- data_bci |> 
#   filter(census > 3)
# 
# data_bci_1_3_unm <- data_unm_fc(data_bci_1_3)
# data_bci_1_3_fil <- data_filter_fc(data_bci_1_3_unm) 
# data_bci_4_8_unm <- data_unm_fc(data_bci_4_8)
# data_bci_4_8_fil <- data_filter_fc(data_bci_4_8_unm) 
# 
# ggplot() + 
#   geom_point(data = data_bci_1_3_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1.5, col="black",inherit.aes = FALSE) 
# 
# data_bci_1_3_fil <- data_bci_1_3_fil |>
#   filter(logDensity >=7)
# 
# mod_lmm_BCI = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|plotID) + (1|species), # + (1|years_since_management_bins),
#                 data = data_bci_1_3_fil, na.action = "na.exclude")
# out <- summary(mod_lmm_BCI)
# out
# plot(allEffects(mod_lmm_BCI))
# plot_model(mod_lmm_BCI,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate=V1, pvalue=V2) |>
#   mutate(estimate = round(estimate,3),
#          pvalue = signif(pvalue,digits=3),
#          pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
# figbci_1_3 <- plot_model(mod_lmm_BCI,type = "pred",show.data=F ,dot.size=1.0,
#                          terms = c("logQMD", "year")) + theme_classic() + 
#   theme(text=element_text(size=10), plot.title=element_text(size=10)) + 
#   geom_point(data = data_bci_1_3_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
#   labs(title = "BCI - census 1-3", 
#        caption = paste0("Year estimate = ", caption$estimate, '\n', "Year p-value = ", caption$pvalue))
# figbci_1_3
# 
# mod_lmm_BCI = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|plotID) + (1|species), # + (1|years_since_management_bins),
#                 data = data_bci_4_8_fil, na.action = "na.exclude")
# out <- summary(mod_lmm_BCI)
# out
# plot(allEffects(mod_lmm_BCI))
# plot_model(mod_lmm_BCI,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate=V1, pvalue=V2) |>
#   mutate(estimate = round(estimate,3),
#          pvalue = signif(pvalue,digits=3),
#          pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
# figbci_4_8 <- plot_model(mod_lmm_BCI,type = "pred",show.data=F ,dot.size=1.0,
#                          terms = c("logQMD", "year")) + theme_classic() + 
#   theme(text=element_text(size=10), plot.title=element_text(size=10)) + 
#   geom_point(data = data_bci_4_8_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
#   labs(title = "BCI - census 4-8", 
#        caption = paste0("Year estimate = ", caption$estimate, '\n', "Year p-value = ", caption$pvalue))
# figbci_4_8

# Interaction models -----------------------------------------------------------

## Biome 1: Tropical & Subtropical Moist Broadleaf Forests  ---------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome1 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome1, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome1)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome1 <- plot_stl_bybiome(
  data_fil_biome1, 
  mod_lmm_int_biome1, 
  name = "Tropical Moist Broadleaf Forests", 
  years = c(1985, 2000, 2015)
) + theme(legend.position = "right")

gg_qmd_int_biome1 <- plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, colors = c("#21918c","#fde725", "#440154"),
                       terms = c("year","logQMD[2,3,4]")) + theme_classic() + ggtitle("Tropical Moist Broadleaf Forests")
gg_qmd_int_biome1

## Biome 4: Temperate Broadleaf & Mixed Forests Forest  ------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome4 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome4, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome4)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome4 <- plot_stl_bybiome(
  data_fil_biome4, 
  mod_lmm_int_biome4, 
  name = "Tropical Moist Broadleaf Forests", 
  years = c(1985, 2000, 2015)
)

## Biome 4: Temperate Broadleaf & Mixed Forests   ------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome4 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome4, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome4)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome4 <- plot_stl_bybiome(
  data_fil_biome4, 
  mod_lmm_int_biome4, 
  name = "Tropical Moist Broadleaf Forests", 
  years = c(1985, 2000, 2015)
)

## Biome 5: Temperate Conifer Forests    ---------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome5 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome5, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome5)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome5 <- plot_stl_bybiome(
  data_fil_biome5, 
  mod_lmm_int_biome5, 
  name = "Tropical Moist Broadleaf Forests", 
  years = c(1985, 2000, 2015)
)

## Biome 6: Boreal Forests/Taiga  ----------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome6 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome6, 
  control = lmerControl(optimizer = "bobyqa"),
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome6)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome6 <- plot_stl_bybiome(
  data_fil_biome6, 
  mod_lmm_int_biome6, 
  name = "Tropical Moist Broadleaf Forests", 
  years = c(1985, 2000, 2015)
)

## Biome 12: Mediterranean Forests  --------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome12 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome12, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome12)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome12 <- plot_stl_bybiome(
  data_fil_biome12, 
  mod_lmm_int_biome12, 
  name = "Tropical Moist Broadleaf Forests", 
  years = c(1985, 2000, 2015)
)

# Publication Figure 1 ---------------------------------------------------------

fig1 <- 
  gg_stl_int_biome1 + 
  gg_stl_int_biome4 + 
  gg_stl_int_biome5 + 
  gg_stl_int_biome6 + 
  gg_stl_int_biome12 +
  guide_area() +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "a",tag_suffix = ")")

ggsave(
  filename = here::here("manuscript/figures/fig1int.png"),
  plot = fig1,
  width = 11, 
  height = 7.5
)

# Quantile regression ----------------------------------------------------------
data_unm <- readRDS(here::here("data/data_unm.rds"))

## Biome 1 Tropical & Subtropical Moist Broadleaf Forests ----------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 1)

# identify disturbed plots
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots)

# plot fraction of disturbed over time
df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), car::logit(fdisturbed))) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = "Year",  y = expression(logit(fraction ~ disturbed)))

ggsave(filename = ".pdf", width = 6, height = 4)

# remove disturbed plots
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

data_unm_biome |> 
  ggplot(aes(logQMD, logDensity, color = year)) +
  geom_point() +
  scale_color_viridis_c()

# no scaling on predictors
fit_lqmm <- lqmm(logDensity ~ logQMD + year,
                 random = ~1,
                 group = plotID,
                 tau = c(0.70, 0.90),
                 data = data_unm_biome,
                 type = "normal"
)

plot_lqmm_bybiome(data_unm_biome, fit_lqmm, name = "Tropical & Subtropical Moist Broadleaf Forests")

### Quantile regression within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)

df_lqmm_byqmdbin |> 
  ggplot(aes(coef_year)) +
  geom_histogram(fill = "grey", color = "black", bins = 12) +
  geom_vline(xintercept = 0.0, linetype = "dotted") +
  theme_classic()

df_lqmm_byqmdbin |> 
  ggplot(aes(bin_lqmm, coef_year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = coef_year_lower, ymax = coef_year_upper), width = 0) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted")

## Biome 4 Temperate Broadleaf & Mixed Forests ---------------------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 4)

# identify disturbed plots
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots)

# plot fraction of disturbed over time
df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), car::logit(fdisturbed))) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = "Year",  y = expression(logit(fraction ~ disturbed)))

# remove disturbed plots
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# no scaling on predictors
fit_lqmm <- lqmm(logDensity ~ logQMD + year,
                 random = ~1,
                 group = plotID,
                 tau = c(0.70, 0.90),
                 data = data_unm_biome,
                 type = "normal"
)

plot_lqmm_bybiome(data_unm_biome, fit_lqmm, name = "Temperate Broadleaf & Mixed Forests")

### Quantile regression within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)

df_lqmm_byqmdbin |> 
  ggplot(aes(coef_year)) +
  geom_histogram(fill = "grey", color = "black", bins = 12) +
  geom_vline(xintercept = 0.0, linetype = "dotted") +
  theme_classic()

df_lqmm_byqmdbin |> 
  ggplot(aes(bin_lqmm, coef_year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = coef_year_lower, ymax = coef_year_upper), width = 0) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted")


## Biome 5  Temperate Conifer Forests Forest ----------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 5)

# identify disturbed plots
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots)

# plot fraction of disturbed over time
df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), car::logit(fdisturbed))) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = "Year",  y = expression(logit(fraction ~ disturbed)))

# remove disturbed plots
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# no scaling on predictors
fit_lqmm <- lqmm(logDensity ~ logQMD + year,
                 random = ~1,
                 group = plotID,
                 tau = c(0.70, 0.90),
                 data = data_unm_biome,
                 type = "normal"
)

plot_lqmm_bybiome(data_unm_biome, fit_lqmm, name = "Temperate Conifer Forests Forest")

### Quantile regression within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)

df_lqmm_byqmdbin |> 
  ggplot(aes(coef_year)) +
  geom_histogram(fill = "grey", color = "black", bins = 12) +
  geom_vline(xintercept = 0.0, linetype = "dotted") +
  theme_classic()

df_lqmm_byqmdbin |> 
  ggplot(aes(bin_lqmm, coef_year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = coef_year_lower, ymax = coef_year_upper), width = 0) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted")


## Biome 6 Boreal Forests/Taiga Forest ----------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 6)

# identify disturbed plots
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots)

# plot fraction of disturbed over time
df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), car::logit(fdisturbed))) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = "Year",  y = expression(logit(fraction ~ disturbed)))

# remove disturbed plots
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# no scaling on predictors
fit_lqmm <- lqmm(logDensity ~ logQMD + year,
                 random = ~1,
                 group = plotID,
                 tau = c(0.70, 0.90),
                 data = data_unm_biome,
                 type = "normal"
)

plot_lqmm_bybiome(data_unm_biome, fit_lqmm, name = "Boreal Forests/Taiga Forest")

### Quantile regression within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)

df_lqmm_byqmdbin |> 
  ggplot(aes(coef_year)) +
  geom_histogram(fill = "grey", color = "black", bins = 12) +
  geom_vline(xintercept = 0.0, linetype = "dotted") +
  theme_classic()

df_lqmm_byqmdbin |> 
  ggplot(aes(bin_lqmm, coef_year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = coef_year_lower, ymax = coef_year_upper), width = 0) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted")


## Biome 12 Mediterranean Forests ----------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 12)

# identify disturbed plots
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots)

# plot fraction of disturbed over time
df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), car::logit(fdisturbed))) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(x = "Year",  y = expression(logit(fraction ~ disturbed)))

# remove disturbed plots
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# no scaling on predictors
fit_lqmm <- lqmm(logDensity ~ logQMD + year,
                 random = ~1,
                 group = plotID,
                 tau = c(0.70, 0.90),
                 data = data_unm_biome,
                 type = "normal"
)

plot_lqmm_bybiome(data_unm_biome, fit_lqmm, name = "Mediterranean Forests")

### Quantile regression within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)

df_lqmm_byqmdbin |> 
  ggplot(aes(coef_year)) +
  geom_histogram(fill = "grey", color = "black", bins = 12) +
  geom_vline(xintercept = 0.0, linetype = "dotted") +
  theme_classic()

df_lqmm_byqmdbin |> 
  ggplot(aes(as.numeric(as.character(bin_lqmm)), coef_year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = coef_year_lower, ymax = coef_year_upper), width = 0) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = expression(log(QMD)))


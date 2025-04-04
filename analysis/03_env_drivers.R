# This script analyses the main environmental drivers affecting the STL changes

# load packages ----
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
library(ingestr)
library(DescTools)
library(corrplot)
library(ggokabeito)

# load data
data_fil_biomes <- readRDS(file.path(here::here(), "/data/inputs/data_fil_biomes.rds"))
data_fil_biomes <- data_fil_biomes |>
  filter(year >= 1980) 

# correlation among variables
M <- as.matrix(data_biomes_fil[,c(27,29,30,33)] %>% distinct())
corrplot(cor(M, use="pairwise.complete.obs"), method="number")

# fit model ----
Fit_Year = lmer(logDensity ~ scale(logQMD) + 
                  scale(year) * scale(ai) + 
                  scale(year) * scale(ndep) + 
                  scale(year) * scale(ORGC) + 
                  scale(year) * scale(PBR) + 
                  (1|dataset/plotID) + (1|species),  
                data = data_fil_biomes, na.action = "na.exclude")
summary(Fit_Year)
r.squaredGLMM(Fit_Year)
AIC(Fit_Year)
plot(allEffects(Fit_Year))
ggplot() + 
  geom_point(data = data_fil_biomes |> filter(country=="Switzerland"),
             aes(x = year, y = ndep, col=dataset), alpha=0.5, size = 1.5, inherit.aes = FALSE) 

plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","ai"))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","ndep"))

plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("year","ai"))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("year","ndep"))

plot_model(Fit_Year)
out <- summary(Fit_Year)

estimates <- round(out$coefficients[,c(1,2,5)],4) |>
  as.data.frame() |>
  rename(pval = `Pr(>|t|)`,
         std = `Std. Error`,
         est = Estimate) |>
  rownames_to_column(var = "var") |>
  mutate(pvalue = ifelse(pval>0.1,"", pval),
         pvalue = ifelse(pval<0.05,"*", pvalue),
         pvalue = ifelse(pval<0.01,"**", pvalue),
         pvalue = ifelse(pval<0.001,"***",pvalue))
estimates
str(estimates)
estimates_plot <- estimates |>
  mutate(eff = ifelse(row_number() %in% 1:7, "Main effect", "Interaction terms"),
         eff = as_factor(eff)) |>
  mutate(varnew = ifelse(var == "scale(year)", "year", var),
         varnew = ifelse(var == "scale(ai)", "MI", varnew),
         varnew = ifelse(var == "scale(ndep)", "Ndep", varnew),
         varnew = ifelse(var == "scale(ORGC)", "ORGC", varnew),
         varnew = ifelse(var == "scale(PBR)", "PBR", varnew),
         varnew = ifelse(var == "scale(year):scale(ai)", "MI", varnew),
         varnew = ifelse(var == "scale(year):scale(ndep)", "Ndep", varnew),
         varnew = ifelse(var == "scale(year):scale(ORGC)", "ORGC", varnew),
         varnew = ifelse(var == "scale(year):scale(PBR)", "PBR", varnew),
         varnew = as_factor(varnew),
         varnew = fct_relevel(varnew,c("PBR","Ndep","ORGC","MI"))) |>
  filter(varnew == "MI"|varnew == "Ndep"|varnew ==  "PBR"|varnew == "ORGC")
str(estimates_plot)

# Figure 2 ----
fig2 <- ggplot(estimates_plot) + 
  geom_bar(aes(y = varnew, weight= est,fill = eff), position = position_stack(reverse = TRUE), width=.5) +
  geom_text(aes(y = varnew,  x= est, label = pvalue), color = "white", size = 5, position = position_stack(vjust = 0.5), vjust = 0.75) + 
  xlab("Coefficients") + ylab("Environmental drivers") + labs(fill = "Predictors") +
  geom_vline(xintercept = 0,color = "black", linetype = "dashed", size = .8) +
  theme_classic() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.text = element_text(size = 10),axis.title = element_text(size = 10),
                           axis.text.y = element_text(hjust = 0.5),
                           legend.text = element_text(size = 9.5),legend.title = element_text(size = 9.5),
                           plot.title = element_text(size = 10),
                           legend.key = element_rect(fill = NA, color = NA),
                           legend.position = c(0.85, 0.15),
                           legend.direction="vertical",
                           legend.box = "horizontal",
                           legend.margin = margin(2, 2, 2, 2),
                           legend.key.size = unit(.6, 'cm'),
                           legend.box.margin = margin(1, 1, 1, 1)) +
  scale_x_continuous(limits = c(-0.04,0.055), breaks = seq(-0.05,0.05,0.025)) +
  scale_y_discrete(labels = c("MI" = "Moisture \nIndex", "ORGC" = "Organic \ncarbon", "Ndep" = "Nitrogen \ndeposition", "PBR" = "Phosporus \navailability")) +
  scale_fill_okabe_ito()
fig2

ggsave(paste0(here::here(), "/manuscript/figures/fig2.png"), width = 8, height = 5, dpi=300)

# This script analyses the changes in the STLs over time (calendar year) by biome.

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
library(here)

# load functions ----
source(here("R/functions.R"))

# load data
data_fil_biomes <- readRDS(here("data/data_fil_biomes.rds"))

plot_map_fil <-  plot_map(data_fil_biomes)
plot_map_fil

plot_stl_fil <- plot_stl(data_fil_biomes)
plot_stl_fil

# Biome 1: Tropical & Subtropical Moist Broadleaf Forests  ----

data_fil_biome1 <- readRDS(here::here("data/data_fil_biome1.rds"))

# # XXX my interpretation because I don't have the file above
# data_fil_biome1 <- data_fil_biomes |> 
#   filter(biomeID == 1)

Fit_Year = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome1, 
  na.action = "na.exclude"
  )

summary(Fit_Year)
r.squaredGLMM(Fit_Year)

plot(allEffects(Fit_Year))
plot_model(Fit_Year)
plot_model(
  Fit_Year,
  type = "pred",
  show.data=TRUE, 
  dot.size=1.0, 
  terms = c("logQMD","year")
  )

out <- summary(Fit_Year)
years <- as.integer(summary(data_fil_biome1$year))

caption <- out$coefficients[3,1] |>
  cbind(out$coefficients[3,5]) |>
  as_tibble() |>
  rename(estimate=V1, pvalue=V2) |>
  mutate(estimate = round(estimate,3),
         pvalue = signif(pvalue,digits=3),
         pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))

plot_model(
  Fit_Year,
  type = "pred",
  show.data=TRUE, 
  dot.size=1.0, 
  colors = c("#21918c","#fde725", "#440154"),
  terms = c("logQMD", "year[1985,2000,2015]")
  ) + 
  theme_classic() 

pred <- ggpredict(
  Fit_Year, 
  terms = c("logQMD","year[1985,2000,2015]"), 
  full.data = TRUE
  ) # full.data = TRUE to include random effects. full.data = FALSE ignores group-specific random effects and gives predictions for an "average" group.

plot(pred, add.data = F) 
preddata <- as.data.frame(pred)

# panel for final plot
fig1_bio1 <- ggplot() + 
  geom_point(data = data_fil_biome1, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
  geom_ribbon(data = preddata, aes(x = x, y = predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=.2,show.legend=T) + 
  geom_smooth(data= preddata, aes(x=x, y=predicted, color=group), method = "lm",fullrange = F,size = .6, se=F) +
  labs(x = "ln QMD", y = "ln N",title = "Tropical Moist Broadleaf Forests",color  = "Year") + #, 
  annotate("text", x = 4.5, y = 9.2, label = paste0("n = ",dim(data_fil_biome1)[1]), 
           size = 4, hjust = 1, vjust = 1, fontface = "italic", color = "gray30") +
       #caption = paste0("n = ",length(na.exclude(data_fil_biome1$year)) ,'\n', "Year estimate = ", caption$estimate, '\n', "Year p-value = ", caption$pvalue)) +
   scale_color_manual("Year", #expression(paste(italic("Year"))), 
                     breaks = c("1985","2000", "2015"), 
                     values = c("#21918c", "#fde725", "#440154")) +
  scale_fill_manual("Year", #expression(paste(italic("Year"))), 
                    breaks = c("1985","2000", "2015"), 
                    values = c("#21918c", "#fde725", "#440154")) +
  theme_classic() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.text = element_text(size = 12),axis.title = element_text(size = 12),
                      legend.text = element_text(size = 10),legend.title = element_text(size = 10),
                      plot.title = element_text(size = 12),
                      legend.key = element_rect(fill = NA, color = NA),
                      legend.position = "bottom",
                      plot.caption = element_text(vjust = -1),
                      plot.title.position = "plot") +
  scale_x_continuous(limits = c(2,4.5),breaks = seq(2,4,1)) +
  scale_y_continuous(limits = c(2.9,9.3),breaks = seq(4,8,2))

fig1_bio1

hist_Year <- ggplot(data_fil_biome1, aes(x=year)) + 
  geom_histogram(color="#FFDB6D", fill="#FFDB6D") + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.margin = unit(c(-.5,.1,.1,.1), "cm")) + 
  ggtitle("") +
  scale_x_continuous("Year", breaks = c(1985,2000,2015)) +
  scale_y_continuous("Frequency", breaks = seq(0,300,50))

hist_Year

## STL shifts ----
# Upward shift
# predict y for a given x
pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1994,1995]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)
preddata <- preddata %>% group_by(x) %>% 
  mutate(upSTL=predicted-lag(predicted)) %>%
  mutate(increment=predicted/lag(predicted)) %>%
  mutate(percent=upSTL*100/lag(predicted))
preddata

change_STL <- preddata %>%
  filter(group=="1995") %>%
  ungroup(x) %>%
  summarise(percent=mean(percent)) %>% 
  pull()
change_STL

# biome 4 ----
# Temperate Broadleaf & Mixed Forests Forest
data_fil_biome4 <- readRDS(file.path(here::here(), "/data/inputs/data_fil_biome4.rds"))

Fit_Year = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins),
  data = data_fil_biome4, 
  na.action = "na.exclude"
  )

summary(Fit_Year)
r.squaredGLMM(Fit_Year)

plot(allEffects(Fit_Year))

plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))

out <- summary(Fit_Year)
years <- as.integer(summary(data_fil_biome4$year))

caption <- out$coefficients[3,1] |>
  cbind(out$coefficients[3,5]) |>
  as_tibble() |>
  rename(estimate=V1, pvalue=V2) |>
  mutate(estimate = round(estimate,3),
         pvalue = signif(pvalue,digits=3),
         pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))

plot_model(
  Fit_Year,
  type = "pred",
  show.data=TRUE,
  dot.size=1.0, 
  colors = c("#21918c","#fde725", "#440154"),
  terms = c("logQMD", "year[1985,2000,2015]")
  ) + 
  theme_classic() 

pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1985,2000,2015]]"), full.data = TRUE)
plot(pred, add.data = F) 

preddata <- as.data.frame(pred)

fig1_bio4 <- ggplot() + 
  geom_point(data = data_fil_biome4, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
  geom_ribbon(data = preddata, aes(x = x, y = predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=.2,show.legend=T) + 
  geom_smooth(data= preddata, aes(x=x, y=predicted, color=group), method = "lm",fullrange = F,size = .6, se=F) +
  labs(x = "ln QMD", y = "ln N",title = "Temperate Broadleaf & Mixed Forests",color  = "Year") + #, 
  annotate("text", x = 4.5, y = 9.2, label = paste0("n = ",dim(data_fil_biome4)[1]), 
           size = 4, hjust = 1, vjust = 1, fontface = "italic", color = "gray30") +
  scale_color_manual("Year", #expression(paste(italic("Year"))), 
                     breaks = c("1985","2000", "2015"), 
                     values = c("#21918c", "#fde725", "#440154")) +
  scale_fill_manual("Year", #expression(paste(italic("Year"))), 
                    breaks = c("1985","2000", "2015"), 
                    values = c("#21918c", "#fde725", "#440154")) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.text = element_text(size = 12),axis.title = element_text(size = 12),
                           legend.text = element_text(size = 10),legend.title = element_text(size = 10),
                           plot.title = element_text(size = 12),
                           legend.key = element_rect(fill = NA, color = NA),
                          legend.position = "bottom",
                          plot.caption = element_text(vjust = 0.1)) +
  scale_x_continuous(limits = c(2,4.5),breaks = seq(2,4,1)) +
  scale_y_continuous(limits = c(2.9,9.3),breaks = seq(4,8,2))
fig1_bio4

## STL shifts ----
# Upward shift
# predict y for a given x
pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1994,1995]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)
preddata <- preddata %>% group_by(x) %>% 
  mutate(upSTL=predicted-lag(predicted)) %>%
  mutate(increment=predicted/lag(predicted)) %>%
  mutate(percent=upSTL*100/lag(predicted))
preddata
change_STL <- preddata %>%
  filter(group=="1995") %>%
  ungroup(x) %>%
  summarise(percent=mean(percent)) %>% pull()
change_STL

# biome 5 ----
# Temperate Conifer Forests Forest
data_fil_biome5 <- readRDS(file.path(here::here(), "/data/inputs/data_fil_biome5.rds"))

Fit_Year = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|dataset/plotID) + (1|species), # + (1|years_since_management_bins),
                data = data_fil_biome5, na.action = "na.exclude")
summary(Fit_Year)
r.squaredGLMM(Fit_Year)
plot(allEffects(Fit_Year))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
out <- summary(Fit_Year)
years <- as.integer(summary(data_fil_biome5$year))
caption <- out$coefficients[3,1] |>
  cbind(out$coefficients[3,5]) |>
  as_tibble() |>
  rename(estimate=V1, pvalue=V2) |>
  mutate(estimate = round(estimate,3),
         pvalue = signif(pvalue,digits=3),
         pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, colors = c("#21918c","#fde725", "#440154"),
                       terms = c("logQMD", "year[1985,2000,2015]")) + theme_classic() 

pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1985,2000,2015]]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)

fig1_bio5 <- ggplot() + 
  geom_point(data = data_fil_biome5, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
  geom_ribbon(data = preddata, aes(x = x, y = predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=.2,show.legend=T) + 
  geom_smooth(data= preddata, aes(x=x, y=predicted, color=group), method = "lm",fullrange = F,size = .6, se=F) +
  labs(x = "ln QMD", y = "ln N",title = "Temperate Conifer Forests",color  = "Year") + #, 
  annotate("text", x = 4.5, y = 9.2, label = paste0("n = ",dim(data_fil_biome5)[1]), 
           size = 4, hjust = 1, vjust = 1, fontface = "italic", color = "gray30") +
  scale_color_manual("Year", #expression(paste(italic("Year"))), 
                     breaks = c("1985","2000", "2015"), 
                     values = c("#21918c", "#fde725", "#440154")) +
  scale_fill_manual("Year", #expression(paste(italic("Year"))), 
                    breaks = c("1985","2000", "2015"), 
                    values = c("#21918c", "#fde725", "#440154")) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          axis.text = element_text(size = 12),axis.title = element_text(size = 12),
                          legend.text = element_text(size = 10),legend.title = element_text(size = 10),
                          plot.title = element_text(size = 12),
                          legend.key = element_rect(fill = NA, color = NA),
                          legend.position = "bottom",
                          plot.caption = element_text(vjust = 0.1)) +
  scale_x_continuous(limits = c(2,4.5),breaks = seq(2,4,1)) +
  scale_y_continuous(limits = c(2.9,9.3),breaks = seq(4,8,2))
fig1_bio5

## STL shifts ----
# Upward shift
# predict y for a given x
pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1994,1995]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)
preddata <- preddata %>% group_by(x) %>% 
  mutate(upSTL=predicted-lag(predicted)) %>%
  mutate(increment=predicted/lag(predicted)) %>%
  mutate(percent=upSTL*100/lag(predicted))
preddata
change_STL <- preddata %>%
  filter(group=="1995") %>%
  ungroup(x) %>%
  summarise(percent=mean(percent)) %>% pull()
change_STL

# biome 6 ----
# Boreal Forests/Taiga Forest
data_fil_biome6 <- readRDS(file.path(here::here(), "/data/inputs/data_fil_biome6.rds"))

Fit_Year = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|dataset/plotID) + (1|species), # + (1|years_since_management_bins),
                data = data_fil_biome6, na.action = "na.exclude")
summary(Fit_Year)
r.squaredGLMM(Fit_Year)
plot(allEffects(Fit_Year))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
out <- summary(Fit_Year)
years <- as.integer(summary(data_fil_biome6$year))
caption <- out$coefficients[3,1] |>
  cbind(out$coefficients[3,5]) |>
  as_tibble() |>
  rename(estimate=V1, pvalue=V2) |>
  mutate(estimate = round(estimate,3),
         pvalue = signif(pvalue,digits=3),
         pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, colors = c("#21918c","#fde725", "#440154"),
                       terms = c("logQMD", "year[1985,2000,2015]")) + theme_classic() 

pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1985,2000,2015]]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)

fig1_bio6 <- ggplot() + 
  geom_point(data = data_fil_biome6, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
  geom_ribbon(data = preddata, aes(x = x, y = predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=.2,show.legend=T) + 
  geom_smooth(data= preddata, aes(x=x, y=predicted, color=group), method = "lm",fullrange = F,size = .6, se=F) +
  labs(x = "ln QMD", y = "ln N",title = "Boreal Forests",color  = "Year") + #, 
  annotate("text", x = 4.5, y = 9.2, label = paste0("n = ",dim(data_fil_biome6)[1]), 
           size = 4, hjust = 1, vjust = 1, fontface = "italic", color = "gray30") +
  scale_color_manual("Year", #expression(paste(italic("Year"))), 
                     breaks = c("1985","2000", "2015"), 
                     values = c("#21918c", "#fde725", "#440154")) +
  scale_fill_manual("Year", #expression(paste(italic("Year"))), 
                    breaks = c("1985","2000", "2015"), 
                    values = c("#21918c", "#fde725", "#440154")) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          axis.text = element_text(size = 12),axis.title = element_text(size = 12),
                          legend.text = element_text(size = 10),legend.title = element_text(size = 10),
                          plot.title = element_text(size = 12),
                          legend.key = element_rect(fill = NA, color = NA),
                          legend.position = "bottom",
                          plot.caption = element_text(vjust = 0.1)) +
  scale_x_continuous(limits = c(2,4.5),breaks = seq(2,4,1)) +
  scale_y_continuous(limits = c(2.9,9.3),breaks = seq(4,8,2))
fig1_bio6

## STL shifts ----
# Upward shift
# predict y for a given x
pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1994,1995]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)
preddata <- preddata %>% group_by(x) %>% 
  mutate(upSTL=predicted-lag(predicted)) %>%
  mutate(increment=predicted/lag(predicted)) %>%
  mutate(percent=upSTL*100/lag(predicted))
preddata
change_STL <- preddata %>%
  filter(group=="1995") %>%
  ungroup(x) %>%
  summarise(percent=mean(percent)) %>% pull()
change_STL

# biome 12 ----
# Mediterranean Forests
data_fil_biome12 <- readRDS(file.path(here::here(), "/data/inputs/data_fil_biome12.rds"))

Fit_Year = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|plotID) + (1|species), # + (1|years_since_management_bins),
                data = data_fil_biome12, na.action = "na.exclude")
summary(Fit_Year)
r.squaredGLMM(Fit_Year)
plot(allEffects(Fit_Year))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
out <- summary(Fit_Year)
years <- as.integer(summary(data_fil_biome12$year))
caption <- out$coefficients[3,1] |>
  cbind(out$coefficients[3,5]) |>
  as_tibble() |>
  rename(estimate=V1, pvalue=V2) |>
  mutate(estimate = round(estimate,3),
         pvalue = signif(pvalue,digits=3),
         pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, colors = c("#21918c","#fde725", "#440154"),
                        terms = c("logQMD", "year[1985,2000,2015]")) + theme_classic()

pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1985,2000,2015]]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)

fig1_bio12 <- ggplot() + 
  geom_point(data = data_fil_biome12, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
  geom_ribbon(data = preddata, aes(x = x, y = predicted,ymin=conf.low,ymax=conf.high,fill=group),alpha=.2,show.legend=T) + 
  geom_smooth(data= preddata, aes(x=x, y=predicted, color=group), method = "lm",fullrange = F,size = .6, se=F) +
  labs(x = "ln QMD", y = "ln N",title = "Mediterranean Forests",color  = "Year") + #, 
  annotate("text", x = 4.5, y = 9.2, label = paste0("n = ",dim(data_fil_biome12)[1]), 
           size = 4, hjust = 1, vjust = 1, fontface = "italic", color = "gray30") +
  scale_color_manual("Year", #expression(paste(italic("Year"))), 
                     breaks = c("1985","2000", "2015"), 
                     values = c("#21918c", "#fde725", "#440154")) +
  scale_fill_manual("Year", #expression(paste(italic("Year"))), 
                    breaks = c("1985","2000", "2015"), 
                    values = c("#21918c", "#fde725", "#440154")) +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          axis.text = element_text(size = 12),axis.title = element_text(size = 12),
                          legend.text = element_text(size = 10),legend.title = element_text(size = 10),
                          plot.title = element_text(size = 12),
                          legend.key = element_rect(fill = NA, color = NA),
                          legend.position = "bottom",
                          plot.caption = element_text(vjust = 0.1)) +
  scale_x_continuous(limits = c(2,4.5),breaks = seq(2,4,1)) +
  scale_y_continuous(limits = c(2.9,9.3),breaks = seq(4,8,2))
fig1_bio12

## STL shifts ----
# Upward shift
# predict y for a given x
pred <- ggpredict(Fit_Year, terms = c("logQMD","year[1994,1995]"), full.data = TRUE)
plot(pred, add.data = F) 
preddata <- as.data.frame(pred)
preddata <- preddata %>% group_by(x) %>% 
  mutate(upSTL=predicted-lag(predicted)) %>%
  mutate(increment=predicted/lag(predicted)) %>%
  mutate(percent=upSTL*100/lag(predicted))
preddata
change_STL <- preddata %>%
  filter(group=="1995") %>%
  ungroup(x) %>%
  summarise(percent=mean(percent)) %>% pull()
change_STL

# Figure 1 ----
fig1 <- fig1_bio1 + fig1_bio4 + fig1_bio5 + fig1_bio6 + fig1_bio12 + guide_area() + plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "a",tag_suffix = ")")
fig1
fig1 <- fig1_bio1 + fig1_bio4 + fig1_bio5 + fig1_bio6 + fig1_bio12 + plot_layout(ncol = 3, guides = "collect") & 
  theme(legend.position = 'bottom') + plot_annotation(tag_levels = "a",tag_suffix = ")")
fig1

ggsave(paste0(here::here(), "/manuscript/figures/fig1.png"), width = 11, height = 7.5, dpi=300)

fig1p <- fig1_bio1 + fig1_bio4 + fig1_bio5 + fig1_bio6 + fig1_bio12  + guide_area() + plot_layout(ncol = 2, guides = "collect") 
fig1p
ggsave(paste0(here::here(), "/manuscript/figures/fig1v.png"), width = 8, height = 12, dpi=300)


# Special case of BCI ----

data_bci <- readRDS(file.path(here::here(), "/data/inputs/data_bci.rds"))
data_bci
plot_stl(data_bci)
plot_map(data_bci)

data_bci_1_3 <- data_bci |> 
  filter(census <= 3)

data_bci_1_3 <- data_bci |> 
  filter(census ==1 |
         census ==2|
         census ==6)


data_bci_4_8 <- data_bci |> 
  filter(census > 3)

data_bci_1_3_unm <- data_unm_fc(data_bci_1_3)
data_bci_1_3_fil <- data_filter_fc(data_bci_1_3_unm) 
data_bci_4_8_unm <- data_unm_fc(data_bci_4_8)
data_bci_4_8_fil <- data_filter_fc(data_bci_4_8_unm) 

ggplot() + 
  geom_point(data = data_bci_1_3_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1.5, col="black",inherit.aes = FALSE) 

data_bci_1_3_fil <- data_bci_1_3_fil |>
  filter(logDensity >=7)

Fit_Year = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|plotID) + (1|species), # + (1|years_since_management_bins),
                data = data_bci_1_3_fil, na.action = "na.exclude")
out <- summary(Fit_Year)
out
plot(allEffects(Fit_Year))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
caption <- out$coefficients[3,1] |>
  cbind(out$coefficients[3,5]) |>
  as_tibble() |>
  rename(estimate=V1, pvalue=V2) |>
  mutate(estimate = round(estimate,3),
         pvalue = signif(pvalue,digits=3),
         pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
figbci_1_3 <- plot_model(Fit_Year,type = "pred",show.data=F ,dot.size=1.0,
                         terms = c("logQMD", "year")) + theme_classic() + 
  theme(text=element_text(size=10), plot.title=element_text(size=10)) + 
  geom_point(data = data_bci_1_3_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
  labs(title = "BCI - census 1-3", 
       caption = paste0("Year estimate = ", caption$estimate, '\n', "Year p-value = ", caption$pvalue))
figbci_1_3

Fit_Year = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|plotID) + (1|species), # + (1|years_since_management_bins),
                data = data_bci_4_8_fil, na.action = "na.exclude")
out <- summary(Fit_Year)
out
plot(allEffects(Fit_Year))
plot_model(Fit_Year,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
caption <- out$coefficients[3,1] |>
  cbind(out$coefficients[3,5]) |>
  as_tibble() |>
  rename(estimate=V1, pvalue=V2) |>
  mutate(estimate = round(estimate,3),
         pvalue = signif(pvalue,digits=3),
         pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
figbci_4_8 <- plot_model(Fit_Year,type = "pred",show.data=F ,dot.size=1.0,
                         terms = c("logQMD", "year")) + theme_classic() + 
  theme(text=element_text(size=10), plot.title=element_text(size=10)) + 
  geom_point(data = data_bci_4_8_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
  labs(title = "BCI - census 4-8", 
       caption = paste0("Year estimate = ", caption$estimate, '\n', "Year p-value = ", caption$pvalue))
figbci_4_8


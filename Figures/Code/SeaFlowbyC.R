library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(fuzzyjoin)
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
library(ggrepel)

# Name your inputs ----
sf_dat_file  <- "MetaData/SeaFlow/SeaFlow_cmap_v1.3_datOnly.csv"
quant_dat_file <- "Intermediates/Quantified_LongDat_Enviro.csv"
quant_cul_fiel <- "Intermediates/Culture_Intermediates/combined_long_withquan.csv"
field_metadat_filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"

# Load up metadat from environmental samples ---
metadat_enviro <- read_csv(field_metadat_filename) %>%
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth, latitude)

# Read homarine dat from environment ----
dat_homarine <- read_csv(quant_dat_file) %>%
  filter(Identification %in% c("Homarine", "Trigonelline")) %>%
  select(Identification, SampID,  nmolCave) %>%
  left_join(metadat_enviro) %>%
  group_by(Station_1, Depth, Cruise, Identification) %>%
  summarise(Lat = mean(latitude),
            nmolCave_mean = mean(nmolCave, na.rm = T), 
            nmolCave_sd = sd(nmolCave, na.rm = T)) %>%
  filter(Cruise == "KOK1606") 
CV.ave <- dat_homarine %>%
  mutate(CV = nmolCave_sd/nmolCave_mean) %>%
  group_by(Identification) %>%
  summarise(CV = mean(CV, na.rm = TRUE))
dat_homarine <- dat_homarine %>%
  left_join(CV.ave) %>%
  mutate(nmolCave_sd = ifelse(is.na(nmolCave_sd), CV*nmolCave_mean, nmolCave_sd))

# Mungde SF data, pull out only KM1606 PC data ---
dat <- read_csv(sf_dat_file) %>%
  filter(cruise == "KOK1606") %>%
  select(time:lon, biomass_synecho, biomass_picoeuk, biomass_prochloro)

dat_summ_syn <- dat %>%
  mutate(lat_bin = round_any(lat , 0.2)) %>%
  mutate(biomass_synecho = biomass_synecho/12) %>% #into umolC per L
  group_by(lat_bin) %>%
  summarise(biomass_synecho_mean = mean(biomass_synecho, na.rm = TRUE),
            biomass_synecho_sd = sd(biomass_synecho, na.rm = TRUE))

dat_summ_pe <- dat %>%
  mutate(lat_bin = round_any(lat , 0.2)) %>%
  mutate(biomass_picoeuk = biomass_picoeuk/12) %>% #into umolC per L
    group_by(lat_bin) %>%
  summarise(biomass_picoeuk_mean = mean(biomass_picoeuk, na.rm = TRUE),
            biomass_picoeuk_sd = sd(biomass_picoeuk, na.rm = TRUE))

dat_summ_prochloro <- dat %>%
  mutate(lat_bin = round_any(lat , 0.2)) %>%
  mutate(biomass_prochloro = biomass_prochloro/12) %>% #into umolC per L
  group_by(lat_bin) %>%
  summarise(biomass_prochloro_mean = mean(biomass_prochloro, na.rm = TRUE),
            biomass_prochloro_sd = sd(biomass_prochloro, na.rm = TRUE))
  
# Plot PC of Syn over space, binned -----
g <- ggplot(data = dat_summ_syn, aes(x = lat_bin, y = biomass_synecho_mean)) +
  geom_ribbon(aes(ymin = biomass_synecho_mean-biomass_synecho_sd, 
                  ymax = biomass_synecho_mean+biomass_synecho_sd))
g  

# Plot PC of picoeuks over space, binned -----
g2 <- ggplot(data = dat_summ_pe, aes(x = lat_bin, y = biomass_picoeuk_mean)) +
  geom_ribbon(aes(ymin = biomass_picoeuk_mean-biomass_picoeuk_sd, 
                  ymax = biomass_picoeuk_mean+biomass_picoeuk_sd))
g2
  
# Plot PC of Pro over space, binned -----
g3 <- ggplot(data = dat_summ_prochloro, aes(x = lat_bin, y = biomass_prochloro_mean)) +
  geom_ribbon(aes(ymin = biomass_prochloro_mean-biomass_prochloro_sd, 
                  ymax = biomass_prochloro_mean+biomass_prochloro_sd))
g3


# Put the best matched syn-PC on homarine and Trigonelline data
dat_homarine_wPC <- difference_left_join(dat_homarine, dat_summ_syn %>%
                                           rename(Lat = lat_bin), max_dist = 0.5,
                                         distance_col = "distance")%>%
  group_by(Lat.x) %>%
  top_n(1, desc(distance)) %>%
  ungroup()

# Try to plot homarine vs syn
g4 <- ggplot(data = dat_homarine_wPC, 
             aes(x = biomass_synecho_mean, y = nmolCave_mean, 
                 label = round_any(Lat.x, .2) )) +
  geom_point()+
  geom_errorbar(aes(ymin = nmolCave_mean - nmolCave_sd, 
                    ymax =  nmolCave_mean + nmolCave_sd))+
  geom_errorbarh(aes(xmin = biomass_synecho_mean - biomass_synecho_sd, 
                    xmax =  biomass_synecho_mean + biomass_synecho_sd))+
  geom_text_repel() +
  facet_wrap(~ Identification, scales = "free") 
  

g4

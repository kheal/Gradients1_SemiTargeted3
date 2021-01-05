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
quant_cul_file <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
field_metadat_filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"

# Get culture dat to get %C as homarine in Syn strains
cul_dat <- read_csv(quant_cul_file) %>%
  filter(Org_Type_Specific == "Synechococcus") %>%
  filter(Identification == "Homarine") %>%
  select(CultureID_short:BioVol_perFilter_uL, Identification, molFractionC_pertotalC) %>%
  group_by(CultureID_short) %>%
  summarise(homarine_molfractionC = mean(molFractionC_pertotalC))

# Load up metadat from environmental samples ---
metadat_enviro <- read_csv(field_metadat_filename) %>%
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth, latitude)

# Read homarine dat from environment ----
dat_homarine <- read_csv(quant_dat_file) %>%
  filter(Identification %in% c("Homarine")) %>%
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
  mutate(biomass_synecho = biomass_synecho/12*1000) %>% #into nmolC per L
  group_by(lat_bin) %>%
  summarise(biomass_synecho_mean = mean(biomass_synecho, na.rm = TRUE),
            biomass_synecho_sd = sd(biomass_synecho, na.rm = TRUE))

# Plot PC of Syn over space, binned -----
g <- ggplot(data = dat_summ_syn, aes(x = lat_bin, y = biomass_synecho_mean)) +
  geom_ribbon(aes(ymin = biomass_synecho_mean-biomass_synecho_sd, 
                  ymax = biomass_synecho_mean+biomass_synecho_sd))
g  

# Put the best matched syn-PC on homarine and Trigonelline data
dat_homarine_wPC <- difference_left_join(dat_homarine, dat_summ_syn %>%
                                           rename(Lat = lat_bin), max_dist = 0.5,
                                         distance_col = "distance")%>%
  group_by(Lat.x) %>%
  top_n(1, desc(distance)) %>%
  ungroup()

# Try to plot homarine vs syn
g4 <- ggplot(data = dat_homarine_wPC, 
             aes(x = biomass_synecho_mean*1000, y = nmolCave_mean, 
                 label = round_any(Lat.x, .2) )) +
  geom_errorbar(aes(ymin = nmolCave_mean - nmolCave_sd, 
                    ymax =  nmolCave_mean + nmolCave_sd), color = "grey")+
  geom_errorbarh(aes(xmin = biomass_synecho_mean*1000 - biomass_synecho_sd*1000, 
                    xmax =  biomass_synecho_mean*1000 + biomass_synecho_sd*1000), height = 0,color = "grey")+
  geom_point()+
#  geom_smooth(method=lm)+
  geom_text_repel() +
# geom_abline(slope = cul_dat$homarine_molfractionC[1], linetype = "dashed")+
#  geom_abline(slope = cul_dat$homarine_molfractionC[2], linetype = "dashed")+
  labs(y= "nmol C homarine L-1",
       x = "nmol C Synechococcus biomass L-1") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "none") 
 # # annotate('text', x = 20000, y = 45, 
 #           label = paste0('homarine as total carbon in Synechoccus 8102 (', 
 #                          round(cul_dat$homarine_molfractionC[1]*100, digits = 2), "%)"), 
 #           angle = 15, size = 3) +
 #  annotate('text', x = 3000, y = 150, 
 #           label = paste0('homarine as \ntotal carbon \nin Synechoccus 7803 \n(', 
 #                          round(cul_dat$homarine_molfractionC[1]*100, digits = 1), "%)"), 
 #           size = 3)


#  facet_wrap(~ Identification, scales = "free") 
g4

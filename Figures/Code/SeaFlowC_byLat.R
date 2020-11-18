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
PC_dat_filename <- "MetaData/PCPN/KOK1606_PCPN_UW_Preliminary_OSU_KRH.csv"
lat_binner <- 0.5

# Load up metadat from environmental samples ---
metadat_enviro <- read_csv(field_metadat_filename) %>%
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth, latitude) 

# Load up PC data
PC_dat <- read_csv(PC_dat_filename) %>%
  mutate(NBorSB = ifelse(row_number() < 132, "NB", "SB"))

PC_dat_summ <- PC_dat %>%
  mutate(lat_bin = round_any(Latitude, 1.5)) %>%
  group_by(lat_bin, NBorSB) %>%
  summarise(biomass_mean = mean(PC, na.rm = TRUE),
            biomass_sd = sd(PC, na.rm = TRUE)) %>%
  mutate(phyto_group = "Total PC")

# Get OG SF dat
#Load Seaflow dat
dat_og <- read_csv(OG_datfile)

# Mungde SF data, pull out only KM1606 PC data ---
# raw SF data biomass is in ug/L, count is cel/uL
dat <- read_csv(sf_dat_file) %>%
  filter(cruise == "KOK1606")%>%
  mutate(NBorSB = ifelse(row_number() < 3164, "NB", "SB")) %>%
  select(time:lon, NBorSB, abundance_prochloro, abundance_synecho, abundance_picoeuk,
         biomass_synecho, biomass_picoeuk, biomass_prochloro) 

dat_summ_syn <- dat %>%
  mutate(lat_bin = round_any(lat , lat_binner)) %>%
  mutate(biomass = biomass_synecho/12) %>% #into umolC per L
  group_by(lat_bin, NBorSB) %>%
  summarise(biomass_mean = mean(biomass, na.rm = TRUE),
            biomass_sd = sd(biomass, na.rm = TRUE),
            abundance_mean = mean(abundance_synecho, na.rm = TRUE),
            abundance_sd = sd(abundance_synecho, na.rm = TRUE)) %>%
  mutate(phyto_group = "Synechococcus")

dat_summ_pe <- dat %>%
  mutate(lat_bin = round_any(lat , lat_binner)) %>%
  mutate(biomass = biomass_picoeuk/12) %>% #into umolC per L
  group_by(lat_bin, NBorSB) %>%
  summarise(biomass_mean = mean(biomass, na.rm = TRUE),
            biomass_sd = sd(biomass, na.rm = TRUE),
            abundance_mean = mean(abundance_picoeuk, na.rm = TRUE),
            abundance_sd = sd(abundance_picoeuk, na.rm = TRUE)) %>%
  mutate(phyto_group = "Picoeukaryotes")

dat_summ_pro <- dat %>%
  mutate(lat_bin = round_any(lat ,lat_binner)) %>%
  mutate(biomass = biomass_prochloro/12) %>% #into umolC per L
  group_by(lat_bin, NBorSB) %>%
  summarise(biomass_mean = mean(biomass, na.rm = TRUE),
            biomass_sd = sd(biomass, na.rm = TRUE),
            abundance_mean = mean(abundance_prochloro, na.rm = TRUE),
            abundance_sd = sd(abundance_prochloro, na.rm = TRUE)) %>%
  mutate(phyto_group = "Prochlorococcus")

# Put the data together
dat_summ_combo <- dat_summ_syn %>%
  rbind(dat_summ_pe) %>%
  rbind(dat_summ_pro) %>%
  rbind(PC_dat_summ) %>%
  filter(NBorSB == "NB") %>%
  mutate(phyto_group = factor(phyto_group, 
                              levels = c("Total PC", "Picoeukaryotes", "Prochlorococcus", "Synechococcus")))


# Plot each group by lat bin
g <- ggplot(data = dat_summ_combo, aes(y = biomass_mean, x = lat_bin))+
  geom_line()+
  geom_ribbon(aes(ymax =  biomass_mean + biomass_sd, 
                    ymin = biomass_mean - biomass_sd), alpha = 0.3) +
  facet_wrap(~  phyto_group, nrow = 4, scales = "free") +
  scale_y_continuous(bquote("Total particulate carbon (\U003BCM L" ^-1*')'),
                     expand = c(0,0)) +
  labs( x = "Latitude") +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        legend.position = "none")
g



save_plot("Figures/Manuscript_figures/Seaflow_lat_withC.pdf", g, base_height = 7, base_width = 6, device=cairo_pdf)

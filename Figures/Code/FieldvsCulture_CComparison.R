library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RColorBrewer)
library(beyonce)
library(RCurl)

#TO DO: check order of magnitude
#Make a double box plot like this: https://stackoverflow.com/questions/46068074/double-box-plots-in-ggplot2

#Name your files -----
quandat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
culture.dat.long.filename <- "Intermediates/Culture_Intermediates/combined_long_withquan.csv"
culture.meta.dat.filename <- "MetaData/CultureMetaData_Ccalcs.csv"
field.meta.dat.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"

#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         BestMatch = Compound.Name_figure) %>%
  select(BestMatch, Identification) %>% unique()

#Load up metadata
meta.dat.enviro <- read_csv(field.meta.dat.filename) %>%
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth)

meta.dat.culture <- read_csv(culture.meta.dat.filename) %>%
  select(CultureID, BioVol_perFilter_uL, nmolC_filtered_final)


#Load up the quan enviro data, get everything in nmol C / C, summarize by each station/depth-----
dat <- read_csv(quandat.file) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)
dat2 <- dat %>%
  select(Identification, SampID,  nmolCave, molFractionC, RankPercent ) %>%
  left_join(meta.dat.enviro) %>%
  filter(Cruise != "MGL1704") %>%
  mutate(nmolmetab_perC = nmolCave*Volume/PC_ave, by = "SampID")
dat3 <- dat2 %>%
  group_by(Station_1, Depth, Cruise, Identification) %>%
  summarise(nmolmetab_perC_enviro = mean(nmolmetab_perC, na.rm = T))
dat4 <- dat3 %>%
  group_by(Identification) %>%
  summarise(nmolmetab_perC_enviro_med = median(nmolmetab_perC_enviro, na.rm = T),
            nmolmetab_perC_enviro_SD = sd(nmolmetab_perC_enviro, na.rm = T)) 



#Load up the culture data, get everythign per C, summarize by each station/depth  ----
dat2cul <- read_csv(culture.dat.long.filename) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)%>%
  select(Identification, MassFeature_Column,ID_rep, PresAbs, intracell_conc_umolCL) %>%
  filter(!is.na(Identification)) %>%
  rename(CultureID  = ID_rep)
dat3cul <- dat2cul  %>%
  left_join(meta.dat.culture, by = "CultureID") %>%
  mutate(nmolmetab_perC = intracell_conc_umolCL*BioVol_perFilter_uL/nmolC_filtered_final, by = "CultureID")
dat4cul <-  dat3cul %>%
  group_by(CultureID, Identification) %>%
  summarise(nmolmetab_perC_cul = mean(nmolmetab_perC, na.rm = T)) %>%
  group_by(Identification) %>%
  summarise(nmolmetab_perC_cul_med = median(nmolmetab_perC_cul, na.rm = T),
            nmolmetab_perC_cul_SD = sd(nmolmetab_perC_cul, na.rm = T)) %>%
  mutate(nmolmetab_perC_cul_med = ifelse(is.na(nmolmetab_perC_cul_med), 10^-4, nmolmetab_perC_cul_med))

#Smash together
dat.combo <- dat4cul %>%
  full_join(dat4, by = "Identification")

#Make a double box plot
g <- ggplot(dat = dat.combo, aes(x = nmolmetab_perC_cul_med, y = nmolmetab_perC_enviro_med, text = Identification)) +
  geom_point() +
  geom_errorbar(aes(ymin = nmolmetab_perC_enviro_med - nmolmetab_perC_enviro_SD,
                    ymax = nmolmetab_perC_enviro_med + nmolmetab_perC_enviro_SD)) + 
  geom_errorbarh(aes(xmin = nmolmetab_perC_cul_med - nmolmetab_perC_cul_SD,
                    xmax = nmolmetab_perC_cul_med + nmolmetab_perC_cul_SD)) + 
  scale_x_log10()+
  scale_y_log10()
g


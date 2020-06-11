library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RColorBrewer)
library(beyonce)
library(RCurl)

#TO DO: make two separate plots - one within ML, one below

#TO DO: get a linear model for the relationship in Log/Log space
#Relatively large uncertainties on independent variables in this study required the use of a maximum likelihood estimate method (York et al., 2004) incorporating bivariate analytical uncertainty for all linear regressions (York,1969; Reed, 1989; York et al., 2004; Cantrell, 2008; Thirumalai et al., 2011), which were performed using published MatlabTM code (Thirumalai et al., 2011). 

#TO DO: separate (by color) by organism type

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
  select(CultureID, BioVol_perFilter_uL, nmolC_filtered_final, Org_Type)

#Load up the quan enviro data, get everything in nmol C / C, summarize by each station/depth-----
dat <- read_csv(quandat.file) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)
dat2 <- dat %>%
  select(Identification, SampID,  nmolCave, molFractionC, RankPercent ) %>%
  left_join(meta.dat.enviro) %>%
  filter(Cruise != "MGL1704") %>%
  mutate(nmolmetab_perC = nmolCave/PC_ave, by = "SampID")
dat3 <- dat2 %>%
  group_by(Station_1, Depth, Cruise, Identification) %>%
  summarise(nmolmetab_perC_enviro = mean(nmolmetab_perC, na.rm = T))
dat4 <- dat3 %>%
  group_by(Identification) %>%
  summarise(nmolmetab_perC_enviro_med = mean(nmolmetab_perC_enviro, na.rm = T),
            nmolmetab_perC_enviro_max = max(nmolmetab_perC_enviro, na.rm = T), 
            nmolmetab_perC_enviro_min = min(nmolmetab_perC_enviro, na.rm = T)) 

#Load up the culture data, get everythin per C, summarize by each station/depth  ----
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
  group_by(CultureID, Identification, Org_Type) %>%
  summarise(nmolmetab_perC_cul = mean(nmolmetab_perC, na.rm = T)) %>%
  group_by(Identification, Org_Type) %>%
  summarise(nmolmetab_perC_cul_med = median(nmolmetab_perC_cul, na.rm = T),
            nmolmetab_perC_cul_max = ifelse(is.na(nmolmetab_perC_cul_med), NA, max(nmolmetab_perC_cul, na.rm = T)), 
            nmolmetab_perC_cul_min = ifelse(is.na(nmolmetab_perC_cul_med), NA, min(nmolmetab_perC_cul, na.rm = T))) %>%
  filter(!is.na(nmolmetab_perC_cul_med)) 
dat5cul <-  dat3cul %>%
  group_by(CultureID, Identification) %>%
  summarise(nmolmetab_perC_cul = mean(nmolmetab_perC, na.rm = T)) %>%
  group_by(Identification) %>%
  summarise(nmolmetab_perC_cul_med = median(nmolmetab_perC_cul, na.rm = T),
            nmolmetab_perC_cul_max = ifelse(is.na(nmolmetab_perC_cul_med), NA, max(nmolmetab_perC_cul, na.rm = T)), 
            nmolmetab_perC_cul_min = ifelse(is.na(nmolmetab_perC_cul_med), NA, min(nmolmetab_perC_cul, na.rm = T))) %>%
  filter(!is.na(nmolmetab_perC_cul_med)) 


#  mutate(nmolmetab_perC_cul_med = ifelse(is.na(nmolmetab_perC_cul_med), 10^-5, nmolmetab_perC_cul_med)) %>%
#  mutate(NotObserved = ifelse(nmolmetab_perC_cul_med == 10^-5, "not observed", NA))

#Smash together
dat.combo <- dat4cul %>%
  full_join(dat4, by = "Identification") %>%
  mutate(nmolmetab_perC_cul_med_Log10 = log10(nmolmetab_perC_cul_med),
         nmolmetab_perC_enviro_med_Log10 = log10(nmolmetab_perC_enviro_med))

dat.combo2 <- dat5cul %>%
  full_join(dat4, by = "Identification") %>%
  mutate(field.cul.fc = nmolmetab_perC_cul_med/nmolmetab_perC_enviro_med)

#Make a plot of each of individual org_types vs ugC, with colors and lines
g <- ggplot(dat = dat.combo, 
          aes(x = nmolmetab_perC_cul_med, y = nmolmetab_perC_enviro_med, color = Org_Type, fill = Org_Type)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE)+
  scale_x_log10()+
  scale_y_log10()

g

#Make a plot of all together
g2 <- ggplot(dat = dat.combo2, 
            aes(x = nmolmetab_perC_cul_med, y = nmolmetab_perC_enviro_med, label =  Identification, 
                color = field.cul.fc < .1 | field.cul.fc > 10)) +
   geom_errorbar(aes(ymin = nmolmetab_perC_enviro_min,
                    ymax = nmolmetab_perC_enviro_max), alpha = 0.3) +
   geom_errorbarh(aes(xmin = nmolmetab_perC_cul_min,
                     xmax = nmolmetab_perC_cul_max), alpha = 0.3) +
  geom_point() +
  geom_abline()+
  scale_x_log10()+
  scale_y_log10()

g2

#save it out
save_plot("Figures/Preliminary/ugCperugC.pdf", g2, base_height = 6, base_width = 10, units = "in")

#for exploration
ggplotly()


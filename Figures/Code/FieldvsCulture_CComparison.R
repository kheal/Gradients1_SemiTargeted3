library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RColorBrewer)
library(beyonce)
library(RCurl)
library(ggrepel)

#TO DO: make two separate plots - one within ML, one below
#TD DO: get rid of Nmar data (not enough surface!)

#Name your files -----
quandat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
culture.dat.long.filename <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
culture.meta.dat.filename <- "MetaData/CultureMetaData.csv"
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
  select(Identification, SampID,  nmolCave, molFractionC, RankPercent, nmolmetab_perC) %>%
  left_join(meta.dat.enviro, by = "SampID") %>%
  filter(Depth < 31)
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
  select(Identification, MassFeature_Column,ID_rep, PresAbs, intracell_conc_umolCL, nmolmetab_perC) %>%
  filter(!is.na(Identification)) %>%
  rename(CultureID  = ID_rep) %>%
  filter(!str_detect(CultureID, "Nmar"))
dat3cul <- dat2cul  %>%
  left_join(meta.dat.culture, by = "CultureID") 
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

#Smash together
dat.combo <- dat4cul %>%
  full_join(dat4, by = "Identification") %>%
  mutate(nmolmetab_perC_cul_med_Log10 = log10(nmolmetab_perC_cul_med),
         nmolmetab_perC_enviro_med_Log10 = log10(nmolmetab_perC_enviro_med))

dat.combo2 <- dat5cul %>%
  full_join(dat4, by = "Identification") %>%
  mutate(Detected = ifelse(is.na(nmolmetab_perC_cul_med), "not detected", "detected")) %>%
  mutate(nmolmetab_perC_cul_med = ifelse(is.na(nmolmetab_perC_cul_med), 4E-5, nmolmetab_perC_cul_med)) %>%
  mutate(field.cul.fc = nmolmetab_perC_cul_med/nmolmetab_perC_enviro_med) 

#Make a plot of all together------
g2 <- ggplot(dat = dat.combo2, 
            aes(x = nmolmetab_perC_cul_med, y = nmolmetab_perC_enviro_med, label =  Identification, 
                color = field.cul.fc < .1 | field.cul.fc > 10)) +
   geom_errorbar(aes(ymin = nmolmetab_perC_enviro_min,
                    ymax = nmolmetab_perC_enviro_max), alpha = 0.3) +
   geom_errorbarh(aes(xmin = nmolmetab_perC_cul_min,
                     xmax = nmolmetab_perC_cul_max), alpha = 0.3) +
  geom_point(aes(shape = Detected), fill = "white") +
  scale_shape_manual(values = c(16, 21))+
  geom_text_repel(dat = dat.combo2 %>% filter(field.cul.fc < .1 | field.cul.fc > 10), 
            aes( x =nmolmetab_perC_cul_med, 
                 y = nmolmetab_perC_enviro_med, label =  Identification ), size = 2.5, fontface = "bold")+
  geom_abline()+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values =c("grey", "black"),)+
  labs(y= "nmol metabolite / umol C (environment)", x = "nmol metabolite / umol C (phytoplankton)") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "none")+
  coord_cartesian(ylim=c(4E-5, 2E1), xlim=c(3E-5, 1E2))
g2

#save it out
save_plot("Figures/Manuscript_figures/ugCperugC.pdf", g2, base_height = 4, base_width = 6, units = "in")

#for exploration
#ggplotly()


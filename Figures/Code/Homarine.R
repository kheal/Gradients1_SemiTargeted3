library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RColorBrewer)
library(beyonce)
library(RCurl)

#TO DO: Make a plot of homarine and tri with structures, latitudinal patterns, organism patterns

#Name your compounds of interest
cmpds <- c("Arsenobetaine", "Glycine betaine")
transform_factor = .01

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
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth, latitude)

meta.dat.culture <- read_csv(culture.meta.dat.filename) %>%
  select(CultureID, BioVol_perFilter_uL, nmolC_filtered_final, Org_Type)

#Load up the quan enviro data-----
dat <- read_csv(quandat.file) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)
dat2 <- dat %>%
  select(Identification, SampID,  nmolinEnviroave, nmolCave, molFractionC, RankPercent ) %>%
  left_join(meta.dat.enviro) %>%
  group_by(Station_1, Depth, Cruise, Identification) %>%
  summarise(Lat = mean(latitude),
            nmolinEnviroave_mean = mean(nmolinEnviroave, na.rm = T), 
            nmolinEnviroave_sd = sd(nmolinEnviroave, na.rm = T)) %>%
  filter(Identification %in% cmpds)

#Plot up depth profiles; latitudinal gradients
g.dpnorth <- ggplot(data = dat2 %>% filter(Cruise == "MGL1704") %>%
                      mutate(nmolinEnviroave_mean = 
                               ifelse(Identification == cmpds[1], nmolinEnviroave_mean, nmolinEnviroave_mean*transform_factor)),
                    aes(x = Depth, y = nmolinEnviroave_mean, color = Identification)) +
  geom_point()+
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ . * transform_factor))+
  scale_x_reverse()+
  coord_flip()
g.dpnorth

g.dpsouth <- ggplot(data = dat2 %>% filter(Cruise == "KM1513") %>%
                      mutate(nmolinEnviroave_mean = 
                               ifelse(Identification == cmpds[1], nmolinEnviroave_mean, nmolinEnviroave_mean*transform_factor),
                             nmolinEnviroave_sd = 
                               ifelse(Identification == cmpds[1], nmolinEnviroave_sd, nmolinEnviroave_sd*transform_factor)),
                    aes(x = Depth, y = nmolinEnviroave_mean, color = Identification)) +
  geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = nmolinEnviroave_mean-nmolinEnviroave_sd, 
                  ymax = nmolinEnviroave_mean+nmolinEnviroave_sd, fill = Identification),
              color = NA, alpha = 0.2)+
  scale_y_continuous(sec.axis = sec_axis(~ . * transform_factor))+
  scale_x_reverse()+
  coord_flip()
g.dpsouth

g.dptransect <- ggplot(data = dat2 %>% filter(Cruise == "KOK1606") %>%
                         mutate(nmolinEnviroave_mean = 
                                  ifelse(Identification == cmpds[1], nmolinEnviroave_mean, nmolinEnviroave_mean*transform_factor),
                                nmolinEnviroave_sd = 
                                  ifelse(Identification == cmpds[1], nmolinEnviroave_sd, nmolinEnviroave_sd*transform_factor)),
                    aes(x = Lat, y = nmolinEnviroave_mean, color = Identification)) +
  geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = nmolinEnviroave_mean-nmolinEnviroave_sd, 
                  ymax = nmolinEnviroave_mean+nmolinEnviroave_sd, fill = Identification),
              color = NA, alpha = 0.2)+
  scale_y_continuous(sec.axis = sec_axis(~ . * transform_factor))
g.dptransect
#TO DO: add ribbons to g.dptransect, need latitude, add in average cv for ones with only n = 1

g.combo <- plot_grid(g.dpnorth, g.dpsouth, g.dptransect, ncol = 1)
g.combo





#Load up the culture data, get everythin per C, summarize by each station/depth  ----
datcul <- read_csv(culture.dat.long.filename) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)%>%
  select(Identification, MassFeature_Column,ID_rep, PresAbs, intracell_conc_umolCL) %>%
  filter(!is.na(Identification)) %>%
  rename(CultureID  = ID_rep)

#Get subset of data







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
                color = field.cul.fc < 5^-2 | field.cul.fc > 10)) +
   geom_errorbar(aes(ymin = nmolmetab_perC_enviro_min,
                    ymax = nmolmetab_perC_enviro_max), alpha = 0.3) +
   geom_errorbarh(aes(xmin = nmolmetab_perC_cul_min,
                     xmax = nmolmetab_perC_cul_max), alpha = 0.3) +
  geom_point() +
  geom_abline()+
  scale_x_log10()+
  scale_y_log10()

g2

ggplotly()


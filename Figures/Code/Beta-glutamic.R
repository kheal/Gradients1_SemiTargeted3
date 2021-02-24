library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RColorBrewer)
library(beyonce)
library(RCurl)
library(magick)

#TO DO: Make a plot of beta-glutamic acid and arsenobetaine with structures, latitudinal patterns, AOA contribution

#Name your compounds of interest
cmpds <- c("b-Glutamic acid", "Arsenobetaine")
transform_factor = 50

#Name your files -----
quandat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
field.meta.dat.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"

#Get list of better names-----
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/4a0d24c8da1366c1a623e462339d13989a14d44f/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         BestMatch = Compound.Name_figure) %>%
  select(BestMatch, Identification) %>% unique()

#Load up metadata-----
meta.dat.enviro <- read_csv(field.meta.dat.filename) %>%
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth, latitude)

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
CV.ave <- dat2 %>%
  mutate(CV = nmolinEnviroave_sd/nmolinEnviroave_mean) %>%
  group_by(Identification) %>%
  summarise(CV = mean(CV, na.rm = TRUE))
dat2 <- dat2 %>%
  left_join(CV.ave) %>%
  mutate(nmolinEnviroave_sd = ifelse(is.na(nmolinEnviroave_sd), CV*nmolinEnviroave_mean, nmolinEnviroave_sd))


#Plot up depth profiles; latitudinal gradients----
g.dpnorth <- ggplot(data = dat2 %>% filter(Cruise == "MGL1704") %>%
                      mutate(nmolinEnviroave_mean = 
                               ifelse(Identification == cmpds[1], nmolinEnviroave_mean, nmolinEnviroave_mean*transform_factor)),
                    aes(x = Depth, y = nmolinEnviroave_mean, shape = Identification)) +
  geom_line() +
  geom_point(fill = "white")+
  scale_y_continuous(sec.axis = sec_axis(~ . / transform_factor, name = "nM Arsenobetaine"))+
  scale_x_reverse()+
  scale_shape_manual(values = c(21, 16))+
  coord_flip()+
  labs(y= expression("nM "*beta*"-glutamic acid"), x = "Depth (m)") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "none") 
g.dpnorth

#save_plot("Figures/Manuscript_figures/BetaGlu_and_Arseno.pdf", g.dpnorth, base_height = 4.5, base_width = 3, units = "in")

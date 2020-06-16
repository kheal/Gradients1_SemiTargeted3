library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RColorBrewer)
library(beyonce)
library(RCurl)
library(magick)

#TO DO: Make a plot of homarine and tri with structures, latitudinal patterns, organism patterns

#Name your compounds of interest
cmpds <- c("Homarine", "Trigonelline")
transform_factor = 50

#Name your files -----
quandat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
culture.dat.long.filename <- "Intermediates/Culture_Intermediates/combined_long_withquan.csv"
culture.meta.dat.filename <- "MetaData/CultureMetaData_Ccalcs.csv"
field.meta.dat.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"
chromat.filename <- "RawOutput/Homarine_EICs.csv"

#Get list of better names-----
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         BestMatch = Compound.Name_figure) %>%
  select(BestMatch, Identification) %>% unique()

#Load up metadata-----
meta.dat.enviro <- read_csv(field.meta.dat.filename) %>%
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth, latitude)

meta.dat.culture <- read_csv(culture.meta.dat.filename) %>%
  select(CultureID, BioVol_perFilter_uL, nmolC_filtered_final, Org_Type)

#Load up chromat data -----
chromat.dat <- read_csv(chromat.filename, skip = 1)

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
                    aes(x = Depth, y = nmolinEnviroave_mean, color = Identification)) +
  geom_point()+
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ . / transform_factor, name = "nM Trigonelline"))+
  scale_x_reverse()+
  scale_color_manual(values = c("#1B9E77", "#B7469B"))+
  coord_flip()+
  labs(y= "nM Homarine", x = "Depth (m)") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "none")

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
  scale_y_continuous(sec.axis = sec_axis(~ . / transform_factor, name = "nM Trigonelline"))+
  scale_x_reverse()+
  scale_color_manual(values = c("#1B9E77", "#B7469B"))+
  scale_fill_manual(values = c("#1B9E77", "#B7469B"))+
  coord_flip() +
  labs(y= "nM Homarine", x = "Depth (m)") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "none")

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
  scale_y_continuous(sec.axis = sec_axis(~ ./ transform_factor, name = "nM Trigonelline")) +
  scale_color_manual(values = c("#1B9E77", "#B7469B"))+
  scale_fill_manual(values = c("#1B9E77", "#B7469B"))+
  labs(y= "nM Homarine", x = "Latitude") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = c(0.3, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, 'lines'))

#Make chromats-----
chromat.dat2 <- chromat.dat %>%
  filter(RT < 11) %>%
  mutate(intensity = ifelse(Column == "RP", intensity*8+1E7, intensity)) %>%
  filter(RT < 5 | Column == "HILIC") %>%
  filter(RT > 4 | Column == "RP")

g.chromat <- ggplot(data =chromat.dat2, 
                    aes(y = intensity, x = RT, color = Column)) +
  geom_line() +
  scale_color_manual(values = c("grey30", "grey60"))+
  labs(y= "Intensity", x = "Retention Time (minutes)") +
  annotate("text", x = c(6.3,8.5), y = c(2.2E8,6E7), label = c("homarine","trigonelline"), 
           fontface = "italic", size = 3)+
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = c(0.1, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 7))
g.chromat


#Combine all the plots -----
g.combo <- plot_grid(g.dptransect, g.dpnorth, g.dpsouth, ncol = 3, rel_widths = c(2,1,1), labels = c("A", "B", "C"))

logo_file <- system.file("extdata", "logo.png", package = "cowplot")
molecule_file <- "Figures/Molecules/Homarine.pdf"
molecule_file2 <- "Figures/Molecules/Trigonelline.pdf"

g.combo2 <- ggdraw(g.combo) + 
  draw_image(molecule_file, x = 0.2, y = 1, hjust = 1, vjust = 1, width = 0.13, height = 0.2)+
  draw_image(molecule_file2, x = 0.2, y = 0.8, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
g.combo2

g.combo3 <- plot_grid(g.combo2, g.chromat, ncol = 1, labels = c("", "D"))
g.combo3

save_plot("Figures/Preliminary/Homarine.pdf", g.combo2, base_height = 8, base_width = 6, units = "in")

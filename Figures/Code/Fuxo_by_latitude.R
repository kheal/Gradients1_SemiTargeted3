library(tidyverse)
library(cowplot)
library(here)
library(plotly)
library(Cairo)
library(limSolve)

#Name your inputs
hplc.dat.file <- "MetaData/HPLCData/KOK1606_HPLC_editted.csv"

#Plot PC and PN vs latitude - maybe add this to the map figure?
dat.all <- read_csv(hplc.dat.file, skip = 1) 

dat <- dat.all %>%
  select(latitude, longitude, `Depth (m)`, `19But`:Zea) %>%
  rename(depth = `Depth (m)`,
         Perid = Peri,
         Fucox = Fuco,
         Pras = Prasino,
         Viol = Viola,
         Allox = Allo,
         Lutein = Lut,
         Zeax = Zea)

dat.long <- dat %>%
  gather(key = "pigment", value = "ngperL", `19But`:Zeax)

dat.long.surface <- dat.long %>%
  filter(depth < 20)

dat.long.surface.summary <- dat.long %>%
  filter(depth < 20) %>%
  mutate(latitude_round = round(latitude, digits = 0)) %>%
  group_by(latitude_round, pigment) %>%
  summarise(ngperLave = mean(ngperL),
            ngperLngsd = sd(ngperL))

#Diatom biomarker
fuco <- ggplot(data = dat.long.surface %>% filter(pigment == "Fucox"), 
               aes(x =latitude, y =  ngperL))+
  geom_point(size = 2)+
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))+
  labs(x = "Latitude", y = "ng Fucoxanthin per L")
fuco

save_plot("Figures/Manuscript_figures/Fucoxanthin.pdf", fuco, base_height = 4, base_width = 4, units = "in")

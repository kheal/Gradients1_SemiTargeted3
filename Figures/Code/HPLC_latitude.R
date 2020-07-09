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

#All separate
c <- ggplot(data = dat.long.surface, aes(x =latitude, y =  ngperL, color = pigment))+
  geom_point(size = 2)+
  facet_wrap(~ pigment, scales = "free")

ggsave("Figures/Preliminary/pigments.pdf", c, width = 10.5, height = 6, units = "in")

#All separate summarized
d <- ggplot(data = dat.long.surface.summary, aes(x =latitude_round, y =  ngperLave))+
 # geom_line(size = 2)+
  geom_ribbon(aes(ymax = ngperLave+ngperLngsd, ymin =ngperLave-ngperLngsd, fill = pigment ), alpha = 0.5)+
  facet_wrap(~ pigment, scales = "free")
d
#Plot as pigment/Chla over space
dat.long.surface.chla <- dat.long.surface %>%
  filter(pigment =="Chla") %>% mutate(Chla = ngperL) %>% select(-pigment, -ngperL)

dat.long.surface.perChla <- dat.long.surface %>%
  filter(pigment !="Chla") %>%
  left_join(dat.long.surface.chla) %>%
  mutate(perchla = ngperL/Chla)

c <- ggplot(data = dat.long.surface.perChla, aes(x =latitude, y =  perchla, color = pigment))+
  geom_point(size = 2)+
  facet_wrap(~ pigment, scales = "free")



#Diatom biomarker
fuco <- ggplot(data = dat.long.surface %>% filter(pigment == "Fucox"), aes(x =latitude, y =  ngperL, color = pigment))+
  geom_point(size = 2)

#

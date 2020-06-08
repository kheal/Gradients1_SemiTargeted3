library(tidyverse)
library(here)


dat.filename <- "Intermediates/Quantified_LongDat_Enviro.csv"
Meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
dat2.filename <- "Intermediates/Longdata.csv"

dat <- read_csv(dat.filename)
meta.dat <- read_csv(Meta.dat.file) %>% select(SampID, latitude, Station_1, Cruise, Depth)
dat2 <- read_csv(dat2.filename)

#Homarine
dat.homarine <- dat %>%
  filter(Identification == "Homarine") %>%
  left_join(meta.dat, by = "SampID")
  
dat.homarine.med <- dat.homarine %>%
  group_by(Identification, Station_1, Cruise, Depth) %>%
  summarise(nmolinEnviroave.med = mean(nmolinEnviroave),
            nmolinEnviroave.sd = sd(nmolinEnviroave),
            latitude = mean(latitude)) %>%
  mutate(nmolinEnviroave.sd = ifelse(is.na(nmolinEnviroave.sd), .2*nmolinEnviroave.med, nmolinEnviroave.sd))

g.homarine.KOK <- ggplot() +
  geom_ribbon(data = dat.homarine.med %>% filter(Cruise == "KOK1606"), 
              aes(x = latitude, ymin = nmolinEnviroave.med-nmolinEnviroave.sd,
                  ymax = nmolinEnviroave.med+nmolinEnviroave.sd), alpha = 0.4)+
    geom_point(data = dat.homarine %>% filter(Cruise == "KOK1606"), 
               aes(x = latitude, y = nmolinEnviroave)) +
  xlab("Latitude") +
  ylab("nmol / L")

g.homarine.KOK

#ggsave("Figures/Presentation_figures/homarineKOK.pdf", g.homarine.KOK, width = 8, height = 5, units = "in")


#Glucosylglycerol
dat.glugly <- dat %>%
  filter(Identification == "Glucosylglycerol") %>%
  left_join(meta.dat, by = "SampID")

dat.glugly.med <- dat.glugly %>%
  group_by(Identification, Station_1, Cruise, Depth) %>%
  summarise(nmolinEnviroave.med = mean(nmolinEnviroave),
            nmolinEnviroave.sd = sd(nmolinEnviroave),
            latitude = mean(latitude)) %>%
  mutate(nmolinEnviroave.sd = ifelse(is.na(nmolinEnviroave.sd), .2*nmolinEnviroave.med, nmolinEnviroave.sd))

g.glugly.KOK <- ggplot() +
  geom_ribbon(data = dat.glugly.med %>% filter(Cruise == "KOK1606"), 
              aes(x = latitude, ymin = nmolinEnviroave.med-nmolinEnviroave.sd,
                  ymax = nmolinEnviroave.med+nmolinEnviroave.sd), alpha = 0.4)+
  geom_point(data = dat.glugly %>% filter(Cruise == "KOK1606"), 
             aes(x = latitude, y = nmolinEnviroave)) +
  xlab("Latitude") +
  ylab("nmol / L")

g.glugly.KOK

#ggsave("Figures/Presentation_figures/glyglyKOK.pdf", g.glugly.KOK,  width = 8, height = 5, units = "in")



#Random MFs
#Try these 4: I162.1125R6.83_HILICPos_X_HILICPos
#Try these 4: I162.1125R8.11_HILICPos_X_HILICPos
#Try these 4: I236.1493R8.23_HILICPos_X_HILICPos
#Try these 4: I236.1492R9.91_HILICPos_X_HILICPos

#First plot of random MF with good pattern-----
MF1 <- "I160.0969R7.56_HILICPos_X_HILICPos"

dat.MF1 <- dat2 %>%
filter(MassFeature_Column == MF1) %>%
  left_join(meta.dat, by = "SampID")

dat.MF1.med <- dat.MF1 %>%
  group_by(MassFeature_Column, Station_1, Cruise, Depth) %>%
  summarise(Adjusted_Area_VolNormed.med = mean(Adjusted_Area_VolNormed),
            Adjusted_Area_VolNormed.sd = sd(Adjusted_Area_VolNormed),
            latitude = mean(latitude)) %>%
  mutate(Adjusted_Area_VolNormed.sd = ifelse(is.na(Adjusted_Area_VolNormed.sd), .2*Adjusted_Area_VolNormed.med, Adjusted_Area_VolNormed.sd))

g.MF1.KOK <- ggplot() +
  geom_ribbon(data = dat.MF1.med %>% filter(Cruise == "KOK1606"), 
              aes(x = latitude, ymin = Adjusted_Area_VolNormed.med-Adjusted_Area_VolNormed.sd,
                  ymax = Adjusted_Area_VolNormed.med+Adjusted_Area_VolNormed.sd), alpha = 0.4)+
  geom_point(data = dat.MF1 %>% filter(Cruise == "KOK1606"), 
             aes(x = latitude, y = Adjusted_Area_VolNormed)) +
  xlab("Latitude") +
  ylab("normalized peak area") +
  ggtitle(MF1)

#Second plot of random MF with good pattern-----
MF2 <- "I162.1125R8.11_HILICPos_X_HILICPos"

dat.MF2 <- dat2 %>%
  filter(MassFeature_Column == MF2) %>%
  left_join(meta.dat, by = "SampID")

dat.MF2.med <- dat.MF2 %>%
  group_by(MassFeature_Column, Station_1, Cruise, Depth) %>%
  summarise(Adjusted_Area_VolNormed.med = mean(Adjusted_Area_VolNormed),
            Adjusted_Area_VolNormed.sd = sd(Adjusted_Area_VolNormed),
            latitude = mean(latitude)) %>%
  mutate(Adjusted_Area_VolNormed.sd = ifelse(is.na(Adjusted_Area_VolNormed.sd), .2*Adjusted_Area_VolNormed.med, Adjusted_Area_VolNormed.sd))

g.MF2.KOK <- ggplot() +
  geom_ribbon(data = dat.MF2.med %>% filter(Cruise == "KOK1606"), 
              aes(x = latitude, ymin = Adjusted_Area_VolNormed.med-Adjusted_Area_VolNormed.sd,
                  ymax = Adjusted_Area_VolNormed.med+Adjusted_Area_VolNormed.sd), alpha = 0.4)+
  geom_point(data = dat.MF2 %>% filter(Cruise == "KOK1606"), 
             aes(x = latitude, y = Adjusted_Area_VolNormed)) +
  xlab("Latitude") +
  ylab("normalized peak area") +
  ggtitle(MF2)

#Third plot of random MF with good pattern-----
MF3 <- "I176.0918R7.43_HILICPos_X_HILICPos"

dat.MF3 <- dat2 %>%
  filter(MassFeature_Column == MF3) %>%
  left_join(meta.dat, by = "SampID")

dat.MF3.med <- dat.MF3 %>%
  group_by(MassFeature_Column, Station_1, Cruise, Depth) %>%
  summarise(Adjusted_Area_VolNormed.med = mean(Adjusted_Area_VolNormed),
            Adjusted_Area_VolNormed.sd = sd(Adjusted_Area_VolNormed),
            latitude = mean(latitude)) %>%
  mutate(Adjusted_Area_VolNormed.sd = ifelse(is.na(Adjusted_Area_VolNormed.sd), .2*Adjusted_Area_VolNormed.med, Adjusted_Area_VolNormed.sd))

g.MF3.KOK <- ggplot() +
  geom_ribbon(data = dat.MF3.med %>% filter(Cruise == "KOK1606"), 
              aes(x = latitude, ymin = Adjusted_Area_VolNormed.med-Adjusted_Area_VolNormed.sd,
                  ymax = Adjusted_Area_VolNormed.med+Adjusted_Area_VolNormed.sd), alpha = 0.4)+
  geom_point(data = dat.MF3 %>% filter(Cruise == "KOK1606"), 
             aes(x = latitude, y = Adjusted_Area_VolNormed)) +
  xlab("Latitude") +
  ylab("normalized peak area") +
  ggtitle(MF3)

g.MF3.KOK

#Fourth plot of random MF with good pattern-----
MF4 <- "I308.1605R5.35_HILICPos_X_HILICPos"

dat.MF4 <- dat2 %>%
  filter(MassFeature_Column == MF4) %>%
  left_join(meta.dat, by = "SampID")

dat.MF4.med <- dat.MF4 %>%
  group_by(MassFeature_Column, Station_1, Cruise, Depth) %>%
  summarise(Adjusted_Area_VolNormed.med = mean(Adjusted_Area_VolNormed),
            Adjusted_Area_VolNormed.sd = sd(Adjusted_Area_VolNormed),
            latitude = mean(latitude)) %>%
  mutate(Adjusted_Area_VolNormed.sd = ifelse(is.na(Adjusted_Area_VolNormed.sd), .2*Adjusted_Area_VolNormed.med, Adjusted_Area_VolNormed.sd))

g.MF4.KOK <- ggplot() +
  geom_ribbon(data = dat.MF4.med %>% filter(Cruise == "KOK1606"), 
              aes(x = latitude, ymin = Adjusted_Area_VolNormed.med-Adjusted_Area_VolNormed.sd,
                  ymax = Adjusted_Area_VolNormed.med+Adjusted_Area_VolNormed.sd), alpha = 0.4)+
  geom_point(data = dat.MF4 %>% filter(Cruise == "KOK1606"), 
             aes(x = latitude, y = Adjusted_Area_VolNormed)) +
  xlab("Latitude") +
  ylab("normalized peak area") +
  ggtitle(MF4)
g.MF4.KOK


#Plot them out
g.MFs.KOK <- plot_grid(g.MF1.KOK, g.MF2.KOK, g.MF3.KOK, g.MF4.KOK, ncol = 1)
g.MFs.KOK

ggsave("Figures/Presentation_figures/randomMFs_KOK.pdf",g.MFs.KOK,  height = 15, width = 8, units = "in")

library(zoo)
library(tidyverse)
library(cowplot)

#Set file names
dat.file.MGL <- "MetaData/CTD_dat/MGL/MGL1704_TM_Cast7.cnv"
dat.file.KM <- "MetaData/CTD_dat/KM/ctd.dat"
dat.file.KOK <- "MetaData/CTD_dat/KOK/G1_CTD_KH.csv"
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"

#Get and tidy the data - temperature, salinity, oxygen, chlorophyll, pressure----
dat.MGL <- read_table("MetaData/CTD_dat/MGL/MGL1704_TM_Cast7.cnv", 
                                          col_names = FALSE, skip = 283)
colnames(dat.MGL) <- c("Scan", "Time", "Pressure", "Temperature1", "Temperature2", "Conductivity1", "Conductivity2", 
               "Salinity1", "Salinity2", "O2", "Fluorescence", "Turbidity", "PAR", "flag")
dat.KM <- read.csv(dat.file.KM) 

dat.KOK <- read.csv(dat.file.KOK)

#Tidy MGL data and KM data-----
dat.MGL.clean <- dat.MGL %>%
  rename(pressure = Pressure,
         temp = Temperature1, 
         salinity = Conductivity1,
         o2 = O2,
         RFU = Fluorescence)%>%
  dplyr::select(pressure, temp, salinity, RFU, o2) %>%
  mutate(index = row_number()) %>%
  filter(index > 13955 & index < 25115) %>% 
  dplyr::select(-index) %>%
  as.matrix() %>%
  rollmean(k = 51) %>%
  as.data.frame()

dat.KM.clean <- dat.KM %>%
  rename(pressure = `CTDPRS`,
         temp = `CTDTMP`, 
         salinity = `CTDSAL`,
         RFU = `CHLPIG`)%>%
  select(pressure, temp, salinity, RFU) %>%
  as.data.frame()

#Tidy up the KOK data------
dat.KOK.sub <- dat.KOK%>%
  filter(matlab.date < 736446.9) %>%
  filter(lat > 23 & lat < 38) %>%
  mutate(lat.round = round(lat, digits = 3)) %>%
  group_by(lat.round) %>%
  summarise(sst = mean(sst, na.rm = T),
            sal = mean(sal, na.rm = T))%>%
  as.matrix() %>%
  rollmedian(k = 3) %>%
  as.data.frame()

test <- ggplot(dat.KOK, aes(x = matlab.date, y = lat))+
  geom_line()
test

#Upload and tidy the KM.samp.data and MGL.samp.data-----
samp.dat <- read_csv(meta.dat.file)
samp.dat.KM <- samp.dat %>%
  filter(str_detect(Cruise, "KM"))

samp.dat.MGL <- samp.dat %>%
  filter(str_detect(Cruise, "MGL"))

samp.dat.KOK <- samp.dat %>%
  filter(str_detect(Cruise, "KOK")) %>%
  filter(!str_detect(SampID, "U1")) %>%
  filter(!str_detect(SampID, "U2")) %>%
  filter(!str_detect(SampID, "S2C1"))


#Make latitude vs salinity and temp plot for KOK------
g.KOK.ctd <- ggplot() +
  geom_line(data = dat.KOK.sub, aes (x = lat.round, y = sst), color = "grey60", size = 1.5)+
  geom_line(data = dat.KOK.sub, aes (x = lat.round, y = (sal*4)-120), color = "blue", size = 1) +
  geom_point(data = samp.dat.KOK, aes (x = latitude, y = 23), shape = 18, size =3)+
  scale_y_continuous(expression(bold(paste("Temperature (", degree*C, ")"))), 
                     sec.axis = sec_axis(~ (. +120)/4, name = "Salinity"))+
  annotate("point", x = 30.5, y = max(dat.KOK.sub$sst), fill = "black", shape =25, size = 2.5)+
  annotate("point", x = 30.5, y = min(dat.KOK.sub$sst), fill = "black", shape =24, size = 2.5)+
  theme(axis.title.y.left = element_text(size = 7, colour = "grey60", face ="bold"),
        axis.title.x = element_text(size = 7),
        axis.title.y.right = element_text(size = 7, colour = "blue", face ="bold"),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        aspect.ratio=0.3)+
       # plot.margin = margin(2, 0.2, 2, 0.2, "cm"))+
  scale_x_continuous("Latitude", limits = c(23, 38), expand = c(0, 0)) 
g.KOK.ctd


#Make a plot for MGL depth profile-------
g.MGL.ctd <- ggplot()+
  geom_line(data = dat.MGL.clean, aes(x = pressure, y = (temp-8)*.8), color = "grey60", size = 1.5)+
  geom_line(data = dat.MGL.clean, aes(x = pressure, y = RFU), color = "chartreuse4", size = 1) +
  geom_point(data = samp.dat.MGL, aes (x = Depth, y = 3.8), shape = 18, size =3)+
  scale_y_continuous("Chla (relative fluorescense units)", limits = c(-.1, 4),
                     sec.axis = sec_axis(~ (. /.8)+8, name = expression(bold(paste("Temperature (", degree*C, ")"))),
                                         breaks = c(6, 7 ,8, 9, 10, 11, 12)), expand = c(0, 0))+
  theme(axis.title.x.top = element_text(size = 7, colour = "grey60", face ="bold"),
        axis.title.x.bottom = element_text(size = 7, colour = "chartreuse4", face ="bold"),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        plot.margin = margin(0.2, 2, .2, 2, "cm"))+
  scale_x_reverse("Depth (m)", limits = c(260, 0), expand = c(0, 0)) +
  coord_flip() 

#Make a plot for KM depth profile-------
g.KM.ctd <- ggplot()+
  geom_line(data = dat.KM.clean, aes(x = pressure, y = (temp-15)*.08), color = "grey60", size = 1.5)+
  geom_line(data = dat.KM.clean, aes(x = pressure, y = RFU), color = "chartreuse4", size = 1) +
    geom_point(data = samp.dat.KM, aes (x = Depth, y = 1), shape = 18, size =3)+
  scale_y_continuous("Chla (relative fluorescense units)", limits = c(0, 1),
                     sec.axis = sec_axis(~ (. /.08)+15, name = expression(bold(paste("Temperature (", degree*C, ")"))),
                                         breaks = c(10 ,15, 20, 25, 30)))+
  scale_x_reverse("Depth (m)", limits = c(250, 0), expand = c(0, 0)) +
  theme(axis.title.x.top = element_text(size = 7, colour = "grey60", face ="bold"),
        axis.title.x.bottom = element_text(size = 7, colour = "chartreuse4", face ="bold"),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        plot.margin = margin(0.2, 2, .2, 2, "cm"))+
  coord_flip() 

detach("package:zoo", unload=TRUE)
rm(list=setdiff(ls(), c("g.MGL.ctd", "g.KM.ctd", "g.KOK.ctd")))

library(tidyverse)
library(ncdf4)
library(raster)
library(maps)
library(mapdata)

#Get Chla data and mean it up
Chla2<- raster("MetaData/ChlaData/A20160892016096.L3m_8D_CHL_chlor_a_9km.nc", var="chlor_a", band =1)
Chla3<- raster("MetaData/ChlaData/A20161052016112.L3m_8D_CHL_chlor_a_9km.nc", var="chlor_a", band =1)
Chla4<- raster("MetaData/ChlaData/A20161132016120.L3m_8D_CHL_chlor_a_9km.nc", var="chlor_a", band =1)
Chla5<- raster("MetaData/ChlaData/A20161212016128.L3m_8D_CHL_chlor_a_9km.nc", var="chlor_a", band =1)
ChlaAve <- mean(Chla2, Chla3, Chla4, Chla5,  na.rm=TRUE)

#convert the rater layer to a data.frame and define the column names
ChlaAve.df<-as.data.frame(rasterToPoints(ChlaAve)) 
colnames(ChlaAve.df)<-c("Longitude","Latitude","ChlaAve")
ChlAve.df.sub <- as_tibble(ChlaAve.df)%>%
  filter(Longitude < -154.5 & Longitude > -162) %>%
  filter(Latitude < 43 & Latitude > 17) %>%
  mutate(ChlaAve.new = ifelse(ChlaAve > 0, ChlaAve, 0))

#get data
world.map <- map_data("worldHires")%>% filter(region == "Hawaii")

detach("package:raster", unload=TRUE)
detach("package:ncdf4", unload=TRUE)
detach("package:mapdata", unload=TRUE)
detach("package:maps", unload=TRUE)


#Get your info on the samples
SampInfo <- read_csv("MetaData/SampInfo_wMetaData_withUTC.csv")
Samps <- SampInfo %>% 
  filter(Day_Extract != 0 | is.na(Day_Extract) ) %>%
  dplyr::select(Cruise, Station, latitude, longitude) %>%
  group_by(Station, Cruise) %>%
  mutate(SampleNumber = n(), 
         longitude = mean (longitude), 
         latitude = mean(latitude)) %>%
  ungroup() %>%  dplyr::select(-Station) %>% unique()

#make your plot
KOK.map<-
  ggplot()+
  geom_raster(data=ChlAve.df.sub,
                       aes(x=Longitude,y=Latitude,fill=ChlaAve))+
  scale_fill_gradient2(low="darkslateblue", mid="olivedrab3", high="olivedrab3",
                       midpoint=0.4, limits = c(0.0, 0.6), na.value = "olivedrab3", breaks = c(0.0, 0.3, 0.6))+
  geom_polygon(data=world.map, aes(x=long, y=lat, group=group), colour="black",size=0.3, fill = "grey60") +
  geom_point(data = Samps, aes(x = longitude, y = latitude, shape = Cruise), show.legend = FALSE) +
  geom_segment(aes(x = -159, y = 30, xend = -159, yend = 24), arrow = arrow(length = unit(0.03, "npc"))  )+
  geom_segment(aes(x = -159, y = 32, xend = -159, yend = 37), arrow = arrow(length = unit(0.03, "npc"))  )+
    scale_y_continuous(expand=c(0,0), limits = c(18,42),breaks=c(20, 25, 30, 35, 40))+
  scale_x_continuous(expand=c(0,0), limits = c(-161.5,-154.5), breaks=c(-160, -158, -156))+
  scale_shape_manual(values = c(8, 16, 8))+
  theme(axis.text.x = element_text(size = 6),
        axis.line.y = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_text(size = 6), 
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.background= element_blank(),
        legend.key = element_rect(colour = NA,fill=NA),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width=unit(0.3,"cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.position = "right",
        aspect.ratio=3,
        plot.margin = margin(0, .1, 0, .1, "cm")) +
  labs(fill = c(expression(paste("Chla \n(mg/m"^"3", ")")))) + 
  annotate("text", x = -159.9, y = 28, label = "NPSG", size = 2.3) +
  annotate("text", x = -159.9, y = 35, label = "NPTZ", size = 2.3) +
  annotate("text", x = -157, y = 41, label = "NPTZ \ndepth profile", size = 2.3) +
  annotate("text", x = -156, y = 23.3, label = "NPSG\ndepth profile", size = 2.3) 
KOK.map

rm(list=setdiff(ls(), c("g.mode.KOK", "g.mode.KM", "g.mode.MGL", "g.MGL.ctd", "g.KM.ctd", "g.tile.KM", "g.tile.MGL", "g.tile.KOK", "g.KOK.ctd", "g.cmps.KOK", "g.cmps.MGL", "g.cmps.KM", "KOK.map")))


library(tidyverse)
library(scales)
library(cowplot)
library(here)
library(ggalluvial)

#Name your inputs
KOK.dat.wide.std.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
MGL.dat.wide.std.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.wide.std.file <- "Intermediates/KM_wide_stand_withclusters.csv"
meta.dat.culture.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"

#######

#Read in your files-----
KOK.dat.long <- read_csv(KOK.dat.wide.std.file)
MGL.dat.long <- read_csv(MGL.dat.wide.std.file)
KM.dat.long <- read_csv(KM.dat.wide.std.file)
meta.culture.dat <- read_csv(meta.dat.culture.filename)

#This line is for making horizontal lines in tile plots ----- 
lindat <- tibble(f=as.factor("1"), x=c(-1,1)*Inf)

#Make a tile plot of the KOK clusters (median reps) -----
KOK.dat.rep.number <-  KOK.dat.long %>% 
  left_join(meta.culture.dat %>% select(SampID, Station_1)) %>% select(SampID, Station_1) %>% unique() %>%
  group_by(Station_1) %>% summarise(reps = n()) %>% ungroup()

KOK.dat.long.med <- KOK.dat.long %>% 
  left_join(read_csv(meta.dat.culture.filename) %>% select(SampID, Station_1)) %>%
  group_by(MassFeature_Column, cluster, cluster_letters, Station_1) %>%
  summarize(std_area = mean(std_area), 
            latitude = mean(latitude) ) %>%
  mutate(latitude = as.factor(round(latitude, digits = 1)))  %>% ungroup %>%
  left_join(KOK.dat.rep.number) %>%
  mutate(std_area = std_area/1) %>%
  mutate(cluster_letters = factor(cluster_letters))%>%
  mutate(std_area = ifelse(std_area == 0 , NA, std_area))
  
#Plot for KOK data ------
g.tile.KOK <- ggplot(dat = KOK.dat.long.med %>% 
                       mutate(std_area = ifelse(std_area > .12, .12, std_area)),
                     aes(x = (latitude), y = MassFeature_Column, fill = std_area, colour="")) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  facet_grid(cluster_letters ~ ., scales="free_y", space="free_y")+
  geom_hline(data= lindat, aes(yintercept=x), col="black", size = 0.5) +
  
  annotate("point", x = 4.5, y = 0, fill = "black", shape =24, size = 3)+
  scale_fill_gradient2(low="cornsilk1", mid="brown4", high="brown4",
                       midpoint=.08, limits = c(0, .12), na.value = "grey80", breaks = c(0.0, 0.04, 0.08, 0.12))  +
  scale_colour_manual(values=c("grey80")) +  
  guides(colour=guide_legend("not \nobserved", override.aes=list(fill="grey80")))+
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 6),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 8, face = "italic", angle=0),
        panel.spacing.y=unit(0, "lines"),
        legend.position = "left",
        legend.justification="bottom",
        legend.margin=margin(0,-15,0,2),
        legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        
        plot.margin = margin(0, 0, 0, 0, "cm"))+
  labs(x ="Latitude", y = "Mass Feature", fill = "Standardized \n peak area")
g.tile.KOK


#Make a tile plot of the MGL clusters (median reps) -----
MGL.dat.long.med <- MGL.dat.long %>% 
  left_join(read_csv(meta.dat.culture.filename) %>% select(SampID, Depth)) %>%
  group_by(MassFeature_Column, cluster, cluster_letters, Depth) %>%
  summarize(std_area = mean(std_area)) %>%
  mutate(Depth = as.factor(round(Depth, digits = 1))) %>% ungroup %>%
  mutate(cluster_letters = factor(cluster_letters)) %>%
  mutate(std_area = ifelse(std_area == 0 , NA, std_area))

g.tile.MGL <- ggplot(dat = MGL.dat.long.med %>% 
                       mutate(std_area = ifelse(std_area > 2.5, 2.5, std_area)) %>%
                       mutate(std_area = ifelse(std_area < -2.5, -2.5, std_area)) %>% filter(cluster_letters != "g"), 
                     aes(x = (MassFeature_Column), y = fct_rev(Depth), fill = std_area, colour="")) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  facet_grid(. ~ cluster_letters, scales="free_x", space="free_x")+
  geom_vline(data= lindat, aes(xintercept=x), col="black", size = 0.5) +
  scale_fill_gradient2(low="slateblue4", mid="cornsilk1", high="brown4",
                       midpoint=0, limits = c(0, 1), na.value = "grey80", breaks = c(0.0, 0.5, 1.0))  +
  scale_colour_manual(values=c("grey80")) +  
  guides(colour=guide_legend("not \nobserved", override.aes=list(fill="grey80")))+
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 6), 
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 8, face = "italic"),
        panel.spacing.x=unit(0, "lines"),
        legend.position = "bottom",
        legend.justification="center",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-15,0,0,0),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        plot.margin = margin(0, 0, 0, 0, "cm"))+
  labs(x ="Mass Feature", y = "Depth (m)", fill = "Standardized \npeak area") 
g.tile.MGL

#Make a tile plot of the KM clusters (median reps) -----
KM.dat.rep.number <-  KM.dat.long %>% 
  left_join(read_csv(meta.dat.culture.filename) %>% select(SampID, Depth)) %>% select(SampID, Depth) %>% unique() %>%
  group_by(Depth) %>% summarise(reps = n()) %>% ungroup()

KM.dat.long.med <- KM.dat.long %>% 
  left_join(read_csv(meta.dat.culture.filename) %>% select(SampID, Depth)) %>%
  group_by(MassFeature_Column, cluster, cluster_letters, Depth) %>%
  summarize(std_area = mean(std_area)) %>% left_join(KM.dat.rep.number)%>%
  mutate(std_area = std_area*as.numeric(1)) %>%
  mutate(Depth = as.factor(round(Depth, digits = 1))) %>% ungroup %>%
  mutate(cluster_letters = factor(cluster_letters))%>%
  mutate(std_area = ifelse(std_area == 0 , NA, std_area))

g.tile.KM <- ggplot(dat = KM.dat.long.med %>% 
                       mutate(std_area = ifelse(std_area > 2.5, 2.5, std_area)) %>%
                       mutate(std_area = ifelse(std_area < -2.5, -2.5, std_area)) %>% filter(cluster_letters != "g"), 
                     aes(x = (MassFeature_Column), y = fct_rev(Depth), fill = std_area, colour="")) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  facet_grid(. ~ cluster_letters, scales="free_x", space="free_x")+
  geom_vline(data= lindat, aes(xintercept=x), col="black", size = 0.5) +
  scale_fill_gradient2(low="slateblue4", mid="cornsilk1", high="brown4",
                       midpoint=0, limits = c(0, 0.4), na.value = "grey80", breaks = c(0.0, 0.2, 0.4))  +
  scale_colour_manual(values=c("grey80")) +  
  guides(colour=guide_legend("not \nobserved", override.aes=list(fill="grey80")))+
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 6), 
        strip.background = element_blank(), 
        strip.text.x = element_text(size = 8, face = "italic"),
        panel.spacing.x=unit(0, "lines"),
        legend.position = "bottom",
        legend.justification="center",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-15,0,0,0),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        plot.margin = margin(0, 0, 0, 0, "cm"))+
  labs(x ="Mass Feature", y = "Depth (m)", fill = "Standardized \npeak area") 
g.tile.KM

rm(list=setdiff(ls(), c("g.tile.KM", "g.tile.MGL", "g.tile.KOK", "g.MGL.ctd", "g.KM.ctd", "g.KOK.ctd")))

library(tidyverse)
library(cowplot)
library(here)
options(readr.num_columns = 0)

#File names
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
MF.dat.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"

#Load files
MGL.dat <- read_csv(MGL.dat.file) #Only 280 compounds
KM.dat <- read_csv(KM.dat.file) #307 compounds
KOK.dat <- read_csv(KOK.dat.file) #Only 318 compounds
Meta.dat <- read_csv(Meta.dat.file)
MF.dat <- read_csv(MF.dat.file)

#Make some mode plots-----
#MGL mode plot -----
MGL.dat.ave.mode <- MGL.dat %>%
  group_by(cluster, cluster_letters, Depth) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)),
            stdev_std_area = as.numeric(sd(std_area)),
            max_area = max(std_area),
            min_area = min(std_area)) %>%
  ungroup() %>%
  mutate(cluster_letters = factor(cluster_letters))

MGL.count <- MGL.dat %>% 
  mutate(cluster_letters = factor(cluster_letters)) %>%
  select(cluster_letters, MassFeature_Column) %>% unique() %>% 
  group_by(cluster_letters) %>% summarise(MF.number = n(),
    MF.percent = round(n()/318*100, digits = 0))

MGL.count.id <- MGL.dat %>%  left_join(MF.dat %>% select(MassFeature_Column, Confidence)) %>%
  mutate(cluster_letters = factor(cluster_letters))%>%
  select(cluster_letters, MassFeature_Column, Confidence) %>% filter(Confidence == 1) %>% unique() %>%
  group_by(cluster_letters) %>% summarise(n = n())

MGL.count <- left_join(MGL.count, MGL.count.id) %>%
  mutate(n = ifelse(is.na(n), 0, n))
#make the MGL plot -----
g.mode.MGL <- ggplot()+
  geom_line(data = MGL.dat %>% filter(cluster_letters != "g"), 
            aes(y = std_area, x = Depth, group = MassFeature_Column), size = 0.5, alpha = 0.2)+
  geom_ribbon(data = MGL.dat.ave.mode %>% filter(cluster_letters != "g"), aes(ymin = ifelse(mean_std_area - stdev_std_area > 0, mean_std_area - stdev_std_area, 0),
                                           ymax = ifelse(mean_std_area + stdev_std_area < 1, mean_std_area + stdev_std_area, 1),
                                           x = Depth), alpha = 0.5, fill = "cornflowerblue") +
  geom_text(data = MGL.count%>% filter(cluster_letters != "g"), aes(x = 100, y = 0.5, 
                                  label = cluster_letters), size = 4, fontface = "italic")+
  geom_text(data = MGL.count%>% filter(cluster_letters != "g"), aes(x = 200, y = 0.5, 
                                  label = paste0(MF.number, " (", MF.percent, "%) ", "MFs; \n", n, " IDd")), size = 2)+
  geom_hline(yintercept = 1/length(unique(MGL.dat$SampID)), linetype = "dashed")+
    scale_y_continuous(sec.axis = dup_axis(), limits = c(0, 1.0), expand = c(0,0), breaks = c(0, 0.2, 0.4, 0.6, 0.8))+
  scale_x_reverse() +
  coord_flip() +
  facet_wrap(cluster_letters ~ ., ncol = 3) +
  xlab("Depth (m)") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
  g.mode.MGL
#save_plot("Figures/Preliminary/MGL_modeplots.pdf", g.mode.MGL, base_height = 8, base_width = 6)

#KM mode plot  ---- 
KM.dat.rep.number <-  KM.dat %>% 
  left_join(Meta.dat %>% select(SampID, Depth)) %>% select(SampID, Depth) %>% unique() %>%
  group_by(Depth) %>% summarise(reps = n()) %>% ungroup()

KM.dat.MF.ave <- KM.dat %>%
  group_by(MassFeature_Column, Depth, cluster, cluster_letters) %>%
  summarise(mean_std_area = as.numeric(mean(std_area))) %>%
  ungroup() %>%
  left_join(KM.dat.rep.number) %>%
  mutate(mean_std_area = mean_std_area * 1) %>%
  mutate(cluster_letters = factor(cluster_letters))

KM.dat.ave.mode <- KM.dat %>%
  group_by(cluster,cluster_letters, Depth) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)),
            stdev_std_area = as.numeric(sd(std_area))) %>%
  ungroup() %>%
  left_join(KM.dat.rep.number) %>%
  mutate(mean_std_area = mean_std_area * 1,
         stdev_std_area = stdev_std_area * 1) %>%
  mutate(cluster_letters = factor(cluster_letters))

KM.count <- KM.dat %>% 
  mutate(cluster_letters = factor(cluster_letters)) %>%
  select(cluster_letters, MassFeature_Column) %>% unique() %>% 
  group_by(cluster_letters) %>% summarise(MF.number = n(),
                                          MF.percent = round(n()/318*100, digits = 0))

KM.count.id <- KM.dat %>%  left_join(MF.dat %>% select(MassFeature_Column, Confidence)) %>%
  mutate(cluster_letters = factor(cluster_letters))%>%
  select(cluster_letters, MassFeature_Column, Confidence) %>% filter(Confidence == 1) %>% unique() %>%
  group_by(cluster_letters) %>% summarise(n = n())

KM.count <- left_join(KM.count, KM.count.id) %>%
  mutate(n = ifelse(is.na(n), 0, n))

#plot KM dat-----
g.mode.KM <- ggplot()+
  geom_line(data = KM.dat.MF.ave, aes(y = mean_std_area, x = Depth, group = MassFeature_Column),  size = 0.5, alpha = 0.2)+
  geom_ribbon(data = KM.dat.ave.mode, aes(ymin = ifelse(mean_std_area - stdev_std_area > 0, mean_std_area - stdev_std_area, 0),
                                                                              ymax = ifelse(mean_std_area + stdev_std_area < .4, 
                                                                                            mean_std_area + stdev_std_area, .4),
                                                                              x = Depth), alpha = 0.5, fill = "cornflowerblue") +
  geom_hline(yintercept = 1/length(unique(KM.dat$SampID)), linetype = "dashed")+
  geom_text(data = KM.count%>% filter(cluster_letters != "g"), aes(x = 50, y = 0.25, 
                                                                    label = cluster_letters), size = 4, fontface = "italic")+
  geom_text(data = KM.count%>% filter(cluster_letters != "g"), aes(x = 60, y = 0.25, 
                                                                    label = paste0(MF.number, " (", MF.percent, "%) ", "MFs; \n", n, " IDd")), size = 2)+
  scale_y_continuous(sec.axis = dup_axis(), limits = c(0, .45), expand = c(0,0), breaks = c(0, 0.2, 0.4))+
  scale_x_reverse() +
  coord_flip() +
  facet_wrap(cluster_letters ~ ., ncol = 2) +
  xlab("Depth (m)") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))

g.mode.KM

#save_plot("Figures/Preliminary/KM_modeplots.pdf", g.mode.KM, base_height = 8, base_width = 6)

#KOK mode plot  ---- 
KOK.dat.MF.ave <- KOK.dat %>%
  left_join(Meta.dat %>% select(SampID, Station_1, latitude))%>%
  group_by(cluster, Station_1, MassFeature_Column, cluster_letters) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)), latitude = mean(latitude)) %>%
  ungroup() %>%
  mutate(cluster_letters = factor(cluster_letters))

KOK.dat.ave.mode <- KOK.dat %>%
  left_join(Meta.dat %>% select(SampID, Station_1))%>%
  group_by(cluster, cluster_letters, Station_1) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)),
            stdev_std_area = as.numeric(sd(std_area)),
            max_area = max(std_area),
            min_area = min(std_area),
            latitude = round(mean(latitude), digits = 1)) %>%
  ungroup() %>%
  mutate(cluster_letters = factor(cluster_letters))

KOK.count <- KOK.dat %>% 
  mutate(cluster_letters = factor(cluster_letters)) %>%
  select(cluster_letters, MassFeature_Column) %>% unique() %>% 
  group_by(cluster_letters) %>% summarise(MF.number = n(),
                                          MF.percent = round(n()/318*100, digits = 0))
  
KOK.count.id <- KOK.dat %>%  left_join(MF.dat %>% select(MassFeature_Column, Confidence)) %>%
  mutate(cluster_letters = factor(cluster_letters))%>%
  select(cluster_letters, MassFeature_Column, Confidence) %>% filter(Confidence == 1) %>% unique() %>%
  group_by(cluster_letters) %>% summarise(n = n())

KOK.count <- left_join(KOK.count, KOK.count.id) %>%
  mutate(n = ifelse(is.na(n), 0, n))

g.mode.KOK <- ggplot()+
  geom_line(data = KOK.dat.MF.ave %>% filter(cluster_letters != "h"), 
            aes(y = ifelse(mean_std_area < 0.22, mean_std_area , 0.22), 
                x = latitude, group = MassFeature_Column), size = 0.5, alpha = 0.2)+
  geom_ribbon(data = KOK.dat.ave.mode %>% filter(cluster_letters != "h"), 
              aes(ymin = ifelse(mean_std_area - stdev_std_area > 0, mean_std_area - stdev_std_area, 0),
                                          ymax = ifelse(mean_std_area + stdev_std_area < 0.22, 
                                                        mean_std_area + stdev_std_area, 0.22),
                                          x = latitude), alpha = 0.5, fill = "cornflowerblue") +
  geom_hline(yintercept = 1/length(unique(KOK.dat$SampID)), linetype = "dashed")+
  geom_text(data = KOK.count %>% filter(cluster_letters != "h"), 
            aes(x = 27, y = .15, label = cluster_letters), size = 4, fontface = "italic")+
  geom_text(data = KOK.count %>% filter(cluster_letters != "h"), 
            aes(x = 27, y = .1, label = paste0(MF.number, " (", MF.percent, "%) ", "MFs; \n", n, " IDd")), size = 2)+
  annotate("point", x = 30.5, y = 0, fill = "black", shape =24, size = 3)+
 # annotate("point", x = 30.5, y = 0.22, fill = "black", shape =25, size = 3)+
  scale_y_continuous(sec.axis = dup_axis(), limits = c(0, 0.22), expand = c(0,0), breaks = c(0, 0.1, 0.2))+
  facet_wrap(cluster_letters ~ ., ncol = 1) +
  xlab("Latitude") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
g.mode.KOK

g.mode.KOK2 <- ggplot()+
  geom_line(data = KOK.dat.MF.ave , 
            aes(y = ifelse(mean_std_area < 0.18, mean_std_area , 0.18), 
                x = latitude, group = MassFeature_Column), size = 1, alpha = 0.05)+
  geom_ribbon(data = KOK.dat.ave.mode, 
              aes(ymin = ifelse(mean_std_area - stdev_std_area > 0, mean_std_area - stdev_std_area, 0),
                  ymax = ifelse(mean_std_area + stdev_std_area < 0.18, 
                                mean_std_area + stdev_std_area, 0.18),
                  x = latitude), alpha = 0.6, fill = "cornflowerblue") +
  geom_hline(yintercept = 1/length(unique(KOK.dat$SampID)), linetype = "dashed")+
  geom_text(data = KOK.count, 
            aes(x = 27, y = .15, label = cluster_letters), size = 4, fontface = "italic")+
  geom_text(data = KOK.count, 
            aes(x = 27, y = .1, label = paste0(MF.number, " (", MF.percent, "%) ", "MFs; \n", n, " IDd")), size = 2)+
  annotate("point", x = 30.5, y = 0, fill = "black", shape =24, size = 3)+
  # annotate("point", x = 30.5, y = 0.22, fill = "black", shape =25, size = 3)+
  scale_y_continuous(limits = c(0, 0.18), expand = c(0,0), breaks = c(0,0.05, 0.1, 0.15))+
  facet_wrap(cluster_letters ~ ., ncol = 1) +
  xlab("Latitude") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
g.mode.KOK2
save_plot("Figures/Presentation_figures/KOK_modeplots.pdf", g.mode.KOK2, base_height = 8, base_width = 3)




#save_plot("Figures/Preliminary/KOK_modeplots.pdf", g.mode.KOK, base_height = 8, base_width = 6)
rm(list=setdiff(ls(), c("g.mode.KOK", "g.mode.KM", "g.mode.MGL", "g.MGL.ctd", "g.KM.ctd", "g.tile.KM", "g.tile.MGL", "g.tile.KOK", "g.KOK.ctd")))


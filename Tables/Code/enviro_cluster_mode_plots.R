library(tidyverse)
library(cowplot)

#Goals
##Triple mode plotwioth all three cruises!
##Alluvial plot connecting MGL1706 to KOK1606 to pull out what oes

#File names
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"

#Load files
MGL.dat <- read_csv(MGL.dat.file)
KM.dat <- read_csv(KM.dat.file)
KOK.dat <- read_csv(KOK.dat.file)



#Make some mode plots-----
#MGL mode plot
MGL.dat.ave.mode <- MGL.dat %>%
  group_by(cluster, Depth) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)),
            stdev_std_area = as.numeric(sd(std_area)),
            max_area = max(std_area),
            min_area = min(std_area)) %>%
  ungroup()

MGL.count <- MGL.dat %>% select(cluster, MassFeature_Column) %>% unique()  %>% 
  group_by(cluster) %>% summarise(n = n())

g.mode.MGL <- ggplot(data = MGL.dat.ave.mode, aes(x = Depth, y = mean_std_area))+
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = mean_std_area - stdev_std_area,
                  ymax = mean_std_area + stdev_std_area), alpha = 0.4) +
  geom_ribbon(aes(ymin = mean_std_area - 2*stdev_std_area,
                  ymax = mean_std_area + 2*stdev_std_area), alpha = 0.2) +
  geom_text(data = MGL.count, aes(x = 100, y = 2, label = paste(n, "compounds")))+
  scale_x_reverse() +
  coord_flip() +
  facet_wrap(cluster ~ ., ncol = 1)

#KM mode plot  ---- 
KM.dat.ave.mode <- KM.dat %>%
  group_by(cluster, Depth) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)),
            stdev_std_area = as.numeric(sd(std_area)),
            max_area = max(std_area),
            min_area = min(std_area)) %>%
  ungroup()

KM.count <- KM.dat %>% select(cluster, MassFeature_Column) %>% unique()  %>% 
  group_by(cluster) %>% summarise(n = n())


g.mode.KM <- ggplot(data = KM.dat.ave.mode, aes(x = Depth, y = mean_std_area))+
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = mean_std_area - stdev_std_area,
                  ymax = mean_std_area + stdev_std_area), alpha = 0.4) +
  geom_ribbon(aes(ymin = mean_std_area - 2*stdev_std_area,
                  ymax = mean_std_area + 2*stdev_std_area), alpha = 0.2) +
  geom_text(data = KM.count, aes(x = 100, y = 2, label = paste(n, "compounds")))+
  scale_x_reverse() +
  coord_flip() +
  facet_wrap(cluster ~ ., ncol = 1)


#KOK mode plot  ---- This doesn't work right now - just the n function doesn't work.....
KOK.dat.ave.mode <- KOK.dat %>%
  mutate(latitude_round = round(latitude, digits = 1)) %>%
  group_by(cluster, latitude_round) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)),
            stdev_std_area = as.numeric(sd(std_area)),
            max_area = max(std_area),
            min_area = min(std_area)) %>%
  ungroup()

KOK.count <- KOK.dat %>% select(cluster, MassFeature_Column) %>% unique()  %>% 
  group_by(cluster) %>% summarise(n = n())


g.mode.KOK <- ggplot(data = KOK.dat.ave.mode, aes(x = latitude_round, y = mean_std_area))+
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = mean_std_area - stdev_std_area,
                  ymax = mean_std_area + stdev_std_area), alpha = 0.4) +
  geom_ribbon(aes(ymin = mean_std_area - 2*stdev_std_area,
                  ymax = mean_std_area + 2*stdev_std_area), alpha = 0.2) +
  geom_text(data = KM.count, aes(x = 30, y = 2, label = paste(n, "compounds")))+
  facet_wrap(cluster ~ ., ncol = 1)






#TO DO - futz with the NAs - put KOK in the middle?
dat.combined <- KM.dat %>% select(MassFeature_Column, cluster) %>%
  rename(cluster.KM = cluster) %>% unique() %>%
  left_join(MGL.dat %>% select(MassFeature_Column, cluster) %>%
              rename(cluster.MGL = cluster) %>% unique())%>%
  left_join(KOK.dat %>% select(MassFeature_Column, cluster) %>%
              rename(cluster.KOK = cluster)%>% unique()) %>%
  group_by(cluster.KM, cluster.MGL, cluster.KOK) %>%
  summarise(Freq = n()) %>% ungroup() %>%
  mutate(cluster.MGL = as.character(cluster.MGL), 
         cluster.KM = as.character(cluster.KM), 
         cluster.KOK = as.character(cluster.KOK))                     

is_alluvia_form(clusters.combined, axes = 1:3, silent = TRUE)

g<- ggplot(data = clusters.combined,
           aes(y = Freq, axis1 = cluster.KM, axis2 = cluster.KOK, axis3 = cluster.MGL)) +
  geom_alluvium(aes(fill = cluster.KOK), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", infer.label = TRUE) +
  scale_x_discrete(limits = c("KM1513", "KOK1606", "MGL1704"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") 

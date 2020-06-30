#Helpful websites: https://briatte.github.io/ggnet/ , https://kateto.net/network-visualization , https://cran.r-project.org/web/packages/ggnetwork/vignettes/ggnetwork.html
#https://igraph.discourse.group/t/error-with-igraph-layout/228/6
#might need an older version of igraph

library(tidyverse) 
library(ggplot2) 
library(plotly)
library(igraph)
library(intergraph)
library(ggnetwork)
library(cowplot)
library(here)
library(ghibli)
library(dutchmasters)
library(wesanderson)
library(vroom)
options(readr.num_columns = 0)

#TO DO: don't plot modes that only have a few mass features in it
#TO DO: shapes that correspond to the # of mass features in each group

#File names
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
MF.dat.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
bootstrap.files <- list.files("Intermediates/BootstrapResults/", full.names = TRUE)

#Load files
MGL.dat <- read_csv(MGL.dat.file) #Only 280 compounds
KM.dat <- read_csv(KM.dat.file) #307 compounds
KOK.dat <- read_csv(KOK.dat.file) #Only 318 compounds
Meta.dat <- read_csv(Meta.dat.file)
MF.dat <- read_csv(MF.dat.file)
bootstrap.results <- vroom(bootstrap.files)


#Make some mode plots-----
#munge data for MGL mode plot -----
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
#make the MGL mode plot -----
MGL.cols.needed <- length(unique(MGL.dat.ave.mode$cluster_letters))
pal.MGL <- c(colorRampPalette(dutchmasters$pearl_earring)(MGL.cols.needed))
g.mode.MGL <- ggplot()+
  geom_ribbon(data = MGL.dat.ave.mode, aes(ymin = ifelse(mean_std_area - stdev_std_area > 0, mean_std_area - stdev_std_area, 0),
                                           ymax = ifelse(mean_std_area + stdev_std_area < 1, mean_std_area + stdev_std_area, 1),
                                           x = Depth, fill = cluster_letters)) +
  geom_text(data = MGL.count, aes(x = 100, y = 0.5, 
                                  label = cluster_letters), size = 4, fontface = "italic")+
  scale_y_continuous(sec.axis = dup_axis(), limits = c(0, 1), expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75, 1.0))+
  scale_x_reverse() +
  scale_fill_manual(values = pal.MGL)+
  coord_flip() +
  facet_wrap(cluster_letters ~ ., ncol = 2) +
  xlab("Depth (m)") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        legend.position = "none")
g.mode.MGL


#munge data for KM mode plot  ---- 
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
KM.cols.needed <- length(unique(KM.dat.ave.mode$cluster_letters))
pal.KM <- c(colorRampPalette(wes_palette("Cavalcanti1"))(KM.cols.needed))
g.mode.KM <- ggplot()+
  geom_ribbon(data = KM.dat.ave.mode, 
              aes(ymin = ifelse(mean_std_area - stdev_std_area > 0,
                                mean_std_area - stdev_std_area, 0),
                  ymax = ifelse(mean_std_area + stdev_std_area < .4,
                                mean_std_area + stdev_std_area, .4),
                  x = Depth, fill = cluster_letters)) +
  geom_text(data = KM.count%>% filter(cluster_letters != "g"), aes(x = 50, y = 0.25, 
                                                                    label = cluster_letters), size = 4, fontface = "italic")+
  scale_y_continuous(sec.axis = dup_axis(), limits = c(0, .45), expand = c(0,0), breaks = c(0, 0.2, 0.4))+
  scale_x_reverse() +
  scale_fill_manual(values = pal.KM)+
  coord_flip() +
  facet_wrap(cluster_letters ~ ., ncol = 1) +
  xlab("Depth (m)") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        legend.position = "none")

g.mode.KM


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

#plot KOK dat-----
KOK.cols.needed <- length(unique(KOK.dat.ave.mode$cluster_letters))
pal.KOK <- c(colorRampPalette(brewer.pal(7,"Dark2"))(KOK.cols.needed)[1:KOK.cols.needed])
g.mode.KOK <- ggplot()+
  geom_ribbon(data = KOK.dat.ave.mode, 
              aes(ymin = ifelse(mean_std_area - stdev_std_area > 0, mean_std_area - stdev_std_area, 0),
                                          ymax = ifelse(mean_std_area + stdev_std_area < 0.22, 
                                                        mean_std_area + stdev_std_area, 0.22),
                                          x = latitude, fill = cluster_letters), alpha = 0.8) +
  geom_text(data = KOK.count, 
            aes(x = 27, y = .15, label = cluster_letters), size = 4, fontface = "italic")+

  scale_y_continuous(sec.axis = dup_axis(), limits = c(0, 0.22), expand = c(0,0), breaks = c(0, 0.1, 0.2))+
  facet_wrap(cluster_letters ~ ., ncol = 2) +
  scale_fill_manual(values = pal.KOK)+
  xlab("Latitude") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        legend.position = "none")
g.mode.KOK

modes.combined <- plot_grid(g.mode.KOK, g.mode.MGL, g.mode.KM, ncol = 3, rel_widths = c(1.3, 1, 0.5))
modes.combined










#Try to plot up a network plot-----
#Munge edges and node data to get them into networkable shape
edges.dat <- bootstrap.results %>%
  separate(cluster_overlap, into = c("node1", "node2"), sep = "&") %>%
  mutate(weight = ifelse(is.na(count), 0, count)) %>%
  mutate(pval = as.numeric(pval),
         significant = pval < 0.11) %>%
  select(node1, node2, weight, significant) %>%
  filter(significant == TRUE)

node.dat <- edges.dat %>% select(node1, weight) %>%
  rbind(edges.dat %>% select(node2, weight) %>% rename(node1 = node2)) %>%
  group_by(node1) %>%
  summarise(count = sum(weight)) %>%
  separate(node1, into = c("dataset", "cluster"), sep = "_", remove = FALSE)

#make an igraph object
net <- graph_from_data_frame(d=edges.dat, vertices=node.dat, directed=FALSE) 

#Barely working from here
#set colors of nodes, and shapes of nodes, size of node
#colrs <- c("gray50", "tomato", "gold", "cornflowerblue")
V(net)$color <- colrs[as.factor(V(net)$dataset)]
shaps <- c("square", "circle", "triangle", "diamond")
V(net)$shapes <- shaps[as.factor(V(net)$dataset)]
V(net)$size <- V(net)$count/10

plot(net)

library(ggraph)
ggraph(net, layout="kk") +
  geom_edge_fan(color="gray50") + 
  geom_node_point(aes(color = dataset,
                      shape = dataset), size = 4) +
  geom_node_text(aes(label = cluster), size=3, color="black", repel=T) +
  theme_void()

#Helpful websites: https://briatte.github.io/ggnet/ , https://kateto.net/network-visualization , https://cran.r-project.org/web/packages/ggnetwork/vignettes/ggnetwork.html
#https://igraph.discourse.group/t/error-with-igraph-layout/228/6
#might need an older version of igraph
source("Figures/Code/SampleMap.R")

library(tidyverse) 
library(ggplot2) 
library(plotly)
library(igraph)
library(intergraph)
library(ggnetwork)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(vroom)
library(ggraph)
options(readr.num_columns = 0)

#TO DO: add clouds with description

#File names
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
MF.dat.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
bootstrap.files <- list.files("Intermediates/BootstrapResults/", full.names = TRUE)
PC.dat.file <- "MetaData/PCPN/DepthProfilesPCPN.csv"

#Load files
MGL.dat <- read_csv(MGL.dat.file) 
KM.dat <- read_csv(KM.dat.file) 
KOK.dat <- read_csv(KOK.dat.file) 
Meta.dat <- read_csv(Meta.dat.file)
MF.dat <- read_csv(MF.dat.file)
bootstrap.results <- vroom(bootstrap.files)

#Set colors-----
MGL.color <- "#DCA258"
KM.color <- "#80A0C7"
KOK.color <- "brown4"
Org.color <- "#00295D"


#Make some mode plots-----
#munge data for MGL mode plot -----
MGL.dat.ave.mode <- MGL.dat %>%
  group_by(cluster, cluster_letters, Depth) %>%
  summarise(mean_std_area = as.numeric(mean(std_area, na.rm = TRUE)),
            stdev_std_area = as.numeric(sd(std_area, na.rm = TRUE)),
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

MGL.dat.ave.mode.2 <- MGL.dat.ave.mode %>%
  left_join(MGL.count, by = "cluster_letters") %>%
  filter(MF.percent > 5) 

#add in the PC data-----
PC.dat <- read_csv(PC.dat.file) %>% head(12) %>%
  filter(str_detect(Cruise, "MGL")) 
#grad the max values for the std clusters
MGL.cluster.maxes <- MGL.dat.ave.mode.2 %>%
  filter(Depth < 35) %>%
  select(cluster_letters, mean_std_area) 
PC.dat.MGL <- list()
for (i in 1:length(MGL.cluster.maxes$cluster_letters)){
  PC.dat.MGL[[i]]<- PC.dat %>%
  mutate(cluster_letters = MGL.cluster.maxes$cluster_letters[i]) %>%
  left_join(MGL.cluster.maxes, by = "cluster_letters") %>%
  mutate(PC_toPlot = Pcaverage_Std*mean_std_area)
}
PC.dat.MGL2 <- do.call(rbind, PC.dat.MGL)


#make the MGL mode plot -----
g.mode.MGL <- ggplot()+
  geom_ribbon(data = MGL.dat.ave.mode.2, 
              aes(ymin = ifelse(mean_std_area - stdev_std_area > 0, mean_std_area - stdev_std_area, 0),
                                           ymax = ifelse(mean_std_area + stdev_std_area < 1, 
                                                         mean_std_area + stdev_std_area, 1),
                                           x = Depth), fill = MGL.color) +
  geom_point(data = PC.dat.MGL2, aes (x = Depth, y = PC_toPlot))+
  geom_text(data = MGL.count %>% filter(MF.percent > 5), 
            aes(x = 100, y = 0.5, label = cluster_letters), size = 4, fontface = "italic")+
  geom_text(data = MGL.count %>% filter(MF.percent > 5), 
            aes(x = 150, y = 0.5, 
                label = paste0(MF.number, " (", MF.percent, "%) ", "MFs; \n", n, " IDd")), size = 2.5)+
  scale_y_continuous(sec.axis = dup_axis(), 
                     limits = c(0, 1.1), expand = c(0,0), breaks = c(0,  0.5,  1.0))+
  scale_x_reverse(limits = c(260,0), expand = c(0,0), breaks = c(0, 100, 200)) +
  coord_flip() +
  facet_wrap(cluster_letters ~ ., ncol = 2) +
  xlab("Depth (m)") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x.top  = element_text(size = 8),
        axis.title.x.bottom  = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x.bottom = element_blank(),
        axis.text.x.top = element_text(size = 7),
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
  summarise(mean_std_area = as.numeric(mean(std_area, na.rm = TRUE)),
            stdev_std_area = as.numeric(sd(std_area, na.rm = TRUE))) %>%
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

#add in the PC data-----
PC.dat <- read_csv(PC.dat.file) %>% head(12) %>%
  filter(str_detect(Cruise, "KM")) 
#grad the max values for the std clusters
KM.cluster.maxes <- KM.dat.ave.mode %>%
  filter(Depth < 35) %>%
  select(cluster_letters, mean_std_area) 
PC.dat.KM <- list()
for (i in 1:length(KM.cluster.maxes$cluster_letters)){
  PC.dat.KM[[i]]<- PC.dat %>%
    mutate(cluster_letters = KM.cluster.maxes$cluster_letters[i]) %>%
    left_join(KM.cluster.maxes, by = "cluster_letters") %>%
    mutate(PC_toPlot = Pcaverage_Std*mean_std_area)
}
PC.dat.KM2 <- do.call(rbind, PC.dat.KM)

#plot KM dat-----
KM.cols.needed <- length(unique(KM.dat.ave.mode$cluster_letters))
pal.KM <- c(colorRampPalette(wes_palette("Cavalcanti1"))(KM.cols.needed))
g.mode.KM <- ggplot()+
  geom_ribbon(data = KM.dat.ave.mode, 
              aes(ymin = ifelse(mean_std_area - stdev_std_area > 0,
                                mean_std_area - stdev_std_area, 0),
                  ymax = ifelse(mean_std_area + stdev_std_area < .4,
                                mean_std_area + stdev_std_area, .4),
                  x = Depth), fill =  KM.color) +
  geom_point(data = PC.dat.KM2, aes (x = Depth, y = PC_toPlot))+
  geom_text(data = KM.count, 
            aes(x = 50, y = 0.25, label = cluster_letters), size = 4, fontface = "italic")+
  geom_text(data = KM.count, 
            aes(x = 85, y = 0.25, 
                label = paste0(MF.number, " (", MF.percent, "%) ", "MFs; \n", n, " IDd")), size = 2.5)+
  scale_y_continuous(sec.axis = dup_axis(), limits = c(0, .45), expand = c(0,0), breaks = c(0, 0.2, 0.4))+
  scale_x_reverse(limits = c(126,0), expand = c(0,0), breaks = c(0, 50, 100)) +
  coord_flip() +
  facet_wrap(cluster_letters ~ ., ncol = 1) +
  xlab("Depth (m)") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.title.x.top  = element_text(size = 8),
        axis.title.x.bottom  = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x.bottom = element_blank(),
        axis.text.x.top = element_text(size = 7),
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
  summarise(mean_std_area = as.numeric(mean(std_area, na.rm = TRUE)),
            stdev_std_area = as.numeric(sd(std_area, na.rm = TRUE)),
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
                                          x = latitude), fill = KOK.color) +
  geom_text(data = KOK.count, 
            aes(x = 30, y = .15, label = cluster_letters), size = 4, fontface = "italic")+
  geom_text(data = KOK.count, 
            aes(x = 30, y = 0.07, 
                label = paste0(MF.number, " (", MF.percent, "%) ", "MFs; \n", n, " IDd")), size = 2.5)+
    scale_y_continuous(sec.axis = dup_axis(), limits = c(0, 0.22), expand = c(0,0), breaks = c(0, 0.1, 0.2))+
  facet_wrap(cluster_letters ~ ., ncol = 2) +
  xlab("Latitude") +
  ylab("Standardized peak area")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y.left = element_text(size = 8),
        axis.title.y.right = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.y.left = element_text(size = 7),
        axis.text.y.right = element_blank(),
        axis.text.x = element_text(size = 7),
        legend.position = "none")
g.mode.KOK

modes.combined <- plot_grid(g.mode.KOK, g.mode.MGL, g.mode.KM, ncol = 3, rel_widths = c(1.3, 1, 0.7), labels = c("A", "B", "C"))
modes.combined



#Try to plot up a network plot-----
#Munge edges and node data to get them into networkable shape
edges.dat <- bootstrap.results %>%
  separate(cluster_overlap, into = c("node1", "node2"), sep = "&") %>%
  mutate(weight = ifelse(is.na(count), 0, count)) %>%
  mutate(pval = as.numeric(pval),
         significant = pval < 0.09) %>%
  select(node1, node2, weight, significant) %>%
  filter(significant == TRUE)

node.dat <- edges.dat %>% select(node1, weight) %>%
  rbind(edges.dat %>% select(node2, weight) %>% rename(node1 = node2)) %>%
  group_by(node1) %>%
  summarise(count = sum(weight)) %>%
  separate(node1, into = c("dataset", "cluster"), sep = "_", remove = FALSE) %>%
  mutate(cluster = cluster %>% str_replace("NA", "not \nobserved" ))

big.nodes <- node.dat %>%
  filter(count > 5)

edges.dat.2 <- edges.dat %>%
  filter(node1 %in% big.nodes$node1) %>%
  filter(node2 %in% big.nodes$node1)

node.dat.2 <- node.dat %>%
  filter(node1 %in% big.nodes$node1) 

#make an igraph object
net <- graph_from_data_frame(d=edges.dat.2, vertices=node.dat.2, directed=FALSE) 

#plot up the network
set.seed(33)
g.net <- ggraph(net, layout = 'fr') + #gem, dh, graphopt, fr, kk, lgl all look decent
  geom_edge_fan(color="gray60", aes()) + 
  geom_node_point(aes(color = dataset,
                      shape = dataset), size = 3)+
  geom_node_text(aes(label = cluster), size=2, color="black", repel=T, fontface = "italic") +
  scale_color_manual(values = c(KM.color, KOK.color, MGL.color, Org.color), labels = c("NPSG depth profile", "Meridional transect", "NPTZ depth profile", "Organisms"))+
  scale_shape_manual(values = c(15,16,17,18), labels = c("NPSG depth profile", "Meridional transect", "NPTZ depth profile", "Organisms"))+
  theme(legend.position = c(0.1, 0.1), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7))
g.net

cloud.file.1 <- "Figures/Annotations/Cloud1.pdf"
cloud.file.2 <- "Figures/Annotations/Cloud2.pdf"
cloud.file.3 <- "Figures/Annotations/Cloud3.pdf"

g.net.2 <- ggdraw() +
  draw_image(cloud.file.1,  x = 0.36, y = -.20, scale = .15) + #cloud for core metabs
  draw_image(cloud.file.3,  x = -.04, y =0.17, scale = .3) + #cloud for dino metabs
  draw_image(cloud.file.1,  x = -.37, y = 0.27, scale = .35) + #cloud for rare metabs
  draw_plot(g.net)+
  draw_label("Core \nMetabolites", x = 0.65, y = 0.3, hjust = 0, fontface = "bold", size = 8) +
  draw_label("Rare \nMetabolites", x = 0.15, y = 0.6, hjust = 0, fontface = "bold", size = 8) +
  draw_label("Dinoflagellate-\nassociated \nMetabolites", x = 0.5, y = 0.75, hjust = 0, fontface = "bold", size = 8) 
g.net.2


#load the map?
nodes.with.map <- plot_grid(KOK.map, g.net.2, ncol = 2, rel_widths = c(1, 2), labels = c("D", "E"))

modes.with.nodes <- plot_grid(modes.combined, nodes.with.map, ncol = 1, rel_heights = c(1, 0.9), labels = c("", "D"))
modes.with.nodes

save_plot("Figures/Manuscript_figures/Modes_Nodes_Map_withPC.pdf", modes.with.nodes, base_height = 8, base_width = 7, units = "in")


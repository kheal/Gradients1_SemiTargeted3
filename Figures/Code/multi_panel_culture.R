library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(nationalparkcolors)
library(here)
options(readr.num_columns = 0)


#Load dat----
org.dat.wide.std.file <- "Intermediates/organs_wide_stand_withclusters.csv"
meta.dat.culture.filename <- "MetaData/CultureMetaData.csv"
MF.dat.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
Quan.dat.file <- "Intermediates/Quantified_MFSummary_Enviro.csv"

meta.dat <- read_csv(meta.dat.culture.filename) %>%
  rename(SampID = CultureID) 
MF.dat <- read_csv(MF.dat.file)
quan.dat <- read_csv(Quan.dat.file)

#Get good order for IDd compounds
dat.quan <- read_csv(Quan.dat.file) %>%
  select(Identification, nmolCmed) %>% arrange(desc(nmolCmed))
MForder <- dat.quan$Identification

#Get decent order for organisms----
org.orders <- meta.dat %>% select(Org_Type_Specific, CultureID_short, Org_Type) %>% unique %>%
  arrange(Org_Type, Org_Type_Specific) %>% filter(!is.na(Org_Type_Specific))

#Load data, turn 0s back into NAs----
dat <- read_csv(org.dat.wide.std.file)

dat.mean <- dat %>% left_join(meta.dat, by = "SampID")%>%
  group_by(MassFeature_Column, CultureID_short, Org_Type, Org_Type_Specific, cluster, cluster_letters) %>%
  summarise(std_area = mean(std_area)) %>%  
  mutate(std_area = ifelse(std_area==0, NA, std_area)) %>% ungroup() %>%
  mutate(CultureID_short = factor(CultureID_short, levels = org.orders$CultureID_short)) 

dat.zeros <- dat.mean %>%
  filter(is.na(std_area)) %>%
  mutate(std_area = 0)

dat.letters <- dat.mean %>% ungroup %>%
  select(cluster_letters) %>% unique()

org.orders <- org.orders %>% mutate(CultureID_short = factor(CultureID_short, levels = org.orders$CultureID_short)) 

#Get information on the number of compounds and ID compounds in each cluster----
Org.MF.dat <- left_join(dat, MF.dat, by = "MassFeature_Column") %>%
  select(MassFeature_Column, cluster, cluster_letters, Confidence, Identification) %>% unique() 

Org.MF.count <- Org.MF.dat %>%
  group_by(cluster_letters) %>%
  summarise(count = n()) %>% ungroup %>%
  mutate(percent.count = round(count/length(Org.MF.dat$MassFeature_Column)*100, digits = 1))

Org.MF.count.IDs <- Org.MF.dat %>%
  filter(Confidence == 1) %>%
  group_by(cluster_letters) %>%
  summarise(count.id = n()) 

Org.MF.count <- left_join(Org.MF.count, Org.MF.count.IDs, by = "cluster_letters")

Org.MF.dat.quan <- Org.MF.dat %>%
  filter(Confidence == 1) %>%
  mutate(Identification = factor(Identification, levels = MForder))



#Plot up a dot plot------
pal <- c(park_palette("Redwoods", 5))
pal2 <- c("deepskyblue4", pal[5], pal[2:4], pal[1])
g.boxplot <- ggplot(data = dat.mean, aes(x = factor(Org_Type), std_area)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.8, binwidth = 0.02, aes(fill = Org_Type, colour =Org_Type ), alpha = 0.7)+
  geom_jitter(data = dat.zeros, aes(x = factor(Org_Type), y = std_area,  fill = Org_Type, colour =Org_Type ),
              shape=16, position=position_jitter(width = 0.3, height = 0.05), alpha = 0.2)+
  geom_text(data = dat.letters, aes(x = 0.5, y = 1.05, label = cluster_letters), size = 2.5, fontface = "italic") +
  geom_text(data = Org.MF.count, aes(x = 2.5, y = 1.05, label = paste0(count, " (", percent.count, "%) ", "MFs; ", count.id, " IDd")), size = 2)+
  facet_wrap(cluster_letters ~ ., ncol = 3)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  scale_y_continuous(limits = c(0, 1.1), expand = c(0,0))+
  scale_fill_manual(values = pal2)+
  scale_color_manual(values = pal2)+
  xlab("Broad organism group") +
  ylab("Standardized peak area")+
  theme(legend.position = "none",
    strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5, angle = -30, vjust = 1, margin = margin(t = -6, r = 0, b = 0, l = 0)),          
    plot.margin = margin(0, 0, 0, 0, "cm"))
g.boxplot


#Plot up tiles-----
g.tileplot <- ggplot(data = dat.mean, aes(x = factor(CultureID_short), y = MassFeature_Column, fill = std_area, colour = "")) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  facet_grid(cluster_letters ~ ., scales="free_y", space="free_y")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, col="black", size = 1.5)+
  scale_fill_gradient2(low="slateblue4", mid="cornsilk1", high="darkolivegreen",
                       midpoint=0, limits = c(0.0, 1), na.value = "grey80", breaks = c(0.0, 0.5, 1.0))+
  scale_colour_manual(values=c("grey80")) +  
  guides(colour=guide_legend("not \nobserved", override.aes=list(fill="grey80")))+
  theme(axis.line.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          strip.background = element_blank(), 
          strip.text.y = element_text(size = 7, face = "italic", angle=0),
          panel.spacing.y=unit(0, "lines"),
          legend.position = "left",
          legend.justification="bottom",
          legend.margin=margin(0,-15,0,2),
          legend.box.margin=margin(0,0,0,0),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          plot.margin = margin(0, .5, 0, .5, "cm"))
g.tileplot

#Make another tile plot to smash together-----
g.tileplot.key <- ggplot(data = org.orders, aes(x = factor(CultureID_short), y = 1, fill = Org_Type), colour = "black") +
  geom_tile() +
  geom_text(aes(label = CultureID_short), angle = -90, size = 1.5, fontface = "italic")+
  scale_fill_manual(values = pal2)+
  xlab("Organisms (colored by broad classification)")+
  theme(axis.line.y = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = -4, r = 0, b = 0, l = 0),size = 6),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none",
        plot.margin = margin(-.1, .5, .2, .5, "cm"))

#Plot up the names of IDd compounds------
g.cmps.org <- ggplot() +
  geom_text(data = Org.MF.dat.quan, aes(y = fct_rev(Identification), x = 1, label = Identification), size = 1.8) +
  facet_wrap(cluster_letters ~ ., ncol = 9, scales = "free") +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = "black",
                                        size = 0.5, linetype = "solid"),
        strip.placement = "inside",
        strip.text.x = element_text(face = "italic", size = 9),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        panel.spacing.x=unit(0, "lines"))
g.cmps.org




#Grid them all together
g.tiles.combo <- plot_grid(g.tileplot, g.tileplot.key, align = "v", axis = "lr", ncol = 1, rel_heights = c(1, 0.15))
g.org.combo <- plot_grid(g.boxplot, g.tiles.combo, ncol = 1, rel_heights = c(1.2, 1), labels = c("A", "B"))
g.org.combo.2 <- plot_grid(g.org.combo, g.cmps.org, ncol = 1, rel_heights = c(1, 0.3), labels = c("", "C"))

g.org.combo.2
save_plot("Figures/Manuscript_figures/Organism_overview.pdf", g.org.combo.2, base_width = 6.5, base_height = 9)

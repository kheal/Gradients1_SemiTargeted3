#TO DO: fix all the joining by errors
#TO DO: check KOK tile plot, why doesn't it match perfectly to the g.clus?
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(ggalluvial)
library(nationalparkcolors)
#library(patchwork)

options(readr.num_columns = 0)

#Name your inputs----
KOK.dat.wide.std.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
org.dat.wide.std.file <- "Intermediates/organs_wide_stand_withclusters.csv"
#org.dat.dat <- "Intermediates/culture_combined_long.csv"
meta.dat.culture.filename <- "MetaData/CultureMetaData.csv"
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
bootstrap.file <- "Intermediates/BootstrapResults/KOKvsOrg.csv"


#This line is for making horizontal lines in tile plots ----- 
lindat <- tibble(f=as.factor("1"), x=c(-1,1)*Inf)

#######

#Read in your files-----
KOK.dat.long <- read_csv(KOK.dat.wide.std.file)
Org.dat.long <- read_csv(org.dat.wide.std.file) 
meta.culture.dat <- read_csv(meta.dat.culture.filename)
bootstrap.dat <- read_csv(bootstrap.file)

#Create alluvial plot from bootstrap output--------
clusters.combined.KOK.org <- KOK.dat.long %>% select(MassFeature_Column, cluster_letters) %>%
  rename(cluster.KOK = cluster_letters) %>% unique() %>%
  left_join(Org.dat.long %>% select(MassFeature_Column, cluster_letters) %>%
              rename(cluster.org = cluster_letters) %>% unique(), by = "MassFeature_Column") %>%
  group_by(cluster.KOK, cluster.org) %>%
  summarise(Freq = n()) %>% ungroup() %>%
  mutate(cluster.org = ifelse(is.na(cluster.org), "NA", cluster.org)) %>%
  mutate(cluster_overlap = paste0("KOK_", cluster.KOK, "&", "Org_", cluster.org)) %>%
  mutate(cluster_overlap = cluster_overlap %>% str_replace("Org_NA", "NA")) %>%
  left_join(bootstrap.dat, by = "cluster_overlap") %>%
  mutate(cluster.org = ifelse(cluster.org == "NA", "not \nobserved", cluster.org)) %>%
  mutate(psig = ifelse(pval == "0.01" | pval == "0.05" | pval == "0.1", "sig", "not sig")) %>%
  mutate(psig = ifelse(is.na(psig), "not sig", psig))

clusters.combined.KOK.org.2 <- clusters.combined.KOK.org %>%
  mutate(cluster.KOK = paste0("KOK_", cluster.KOK))

clusters.combined.KOK.org.3 <- clusters.combined.KOK.org.2 %>%
  select(Freq, cluster.KOK, cluster.org, psig) %>% as.data.frame()

#Make the alluvial plot------
#at some point, this worked, now its not :(
gclus<- ggplot(data = clusters.combined.KOK.org.2,
               aes(y = Freq, axis1 = cluster.KOK, axis2 = cluster.org)) +
  geom_alluvium(aes(alpha = psig), 
                fill = "black", width = 1/12,  lode.guidance = "forward") +
 #  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  scale_x_discrete(limits = c("KOK1606_clusters", "Organism_clusters"), 
                   expand = c(0, 0)) +
#  scale_y_discrete(limits = c(0, 1)) +
  scale_alpha_manual(values = c("0.1",  "1"))+
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim=c(15,299), xlim=c(0.95,2.05))

gclus

#Tidy the data for tile plots - the mass features in KOK dat----
KOK.MF.order <- KOK.dat.long %>%
  select(MassFeature_Column, cluster_letters) %>%
  left_join(Org.dat.long %>% select(MassFeature_Column, cluster_letters) %>% rename(Org_clusters = cluster_letters)) %>%
  unique() %>% arrange(cluster_letters, Org_clusters)

Org.MF.order <- KOK.dat.long %>%
  select(MassFeature_Column, cluster_letters) %>%
  left_join(Org.dat.long %>% select(MassFeature_Column, cluster_letters) %>% rename(Org_clusters = cluster_letters)) %>%
  unique() %>% arrange(Org_clusters, cluster_letters)



#Make a tile plot of the KOK medians-------
KOK.dat.long.med <- KOK.dat.long %>% 
  left_join(read_csv(meta.dat.file) %>% select(SampID, Station_1)) %>%
  group_by(MassFeature_Column, cluster, Station_1, cluster_letters) %>%
  summarize(std_area = mean(std_area), 
            latitude = mean(latitude) ) %>%
  mutate(std_area = ifelse(std_area == 0, NA, std_area)) %>%
  mutate(latitude = as.numeric(round(latitude, digits = 1))) 

KOK.dat.long.med$MassFeature_Column <- factor(KOK.dat.long.med$MassFeature_Column, 
                                              levels = rev(KOK.MF.order$MassFeature_Column))
KOK.dat.long.med$latitude <- factor(KOK.dat.long.med$latitude, 
                                              levels = sort(unique(KOK.dat.long.med$latitude)))

g.tile.KOK <- ggplot(dat = KOK.dat.long.med %>% 
                       mutate(std_area = ifelse(std_area > .12, .12, std_area)),
                     aes(x = (latitude), y = MassFeature_Column, fill = std_area, colour="")) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  facet_grid(cluster_letters ~ ., scales="free_y", space="free_y", switch = "y")+
  geom_hline(data= lindat, aes(yintercept=x), col="black", size = 0.5) +
  annotate("point", x = 4.5, y = 0, fill = "black", shape =24, size = 3)+
  scale_fill_gradient2(low="cornsilk1", mid="brown4", high="brown4",
                       midpoint=.08, limits = c(0, .12), na.value = "grey80", breaks = c(0.0, 0.04, 0.08, 0.12))  +
  scale_colour_manual(values=c("grey90")) +  
  guides(colour=guide_legend("not \nobserved", override.aes=list(fill="grey90")))+
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 6),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
    #    strip.text.placement = "outside",
        strip.background = element_blank(), 
        strip.text.y.left = element_text(size = 8, face = "italic", angle=0),
        panel.spacing.y=unit(0, "lines"),
        legend.position = "top",
        legend.justification="left",
        legend.margin=margin(0, 0 ,-15,0),
        legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        
        plot.margin = margin(.25, 0, 0, 0, "cm"))+
  labs(x ="Latitude", y = "Mass Feature", fill = "Standardized \n peak area")
g.tile.KOK

test <- KOK.dat.long.med %>% filter(is.na(cluster_letters))


#Tidy data for tile plot (organisms)-----
orgs.dat.long.med <- Org.dat.long %>% 
  left_join(meta.culture.dat %>% rename(SampID = CultureID)) %>%
  group_by(CultureID_short, MassFeature_Column, Org_Type, Org_Type_Specific, cluster, cluster_letters) %>%
  summarise(std_area = mean(std_area)) %>%
  mutate(std_area = ifelse(std_area == 0 , NA, std_area))

#Add in MFs that we didn't see in oganisms
mfs.with.NA.in.orgs <- KOK.dat.long %>%
  select(MassFeature_Column, cluster_letters)%>%
  rename(cluster.KOK = cluster_letters) %>% unique() %>%
  left_join(orgs.dat.long.med %>% ungroup() %>% select(MassFeature_Column, cluster_letters) %>%
              rename(cluster.org = cluster_letters), by = "MassFeature_Column") %>%
  filter(is.na(cluster.org)) %>% unique() %>%
  mutate(Match = "match")

matcher <- orgs.dat.long.med %>% ungroup() %>%
  select(CultureID_short, Org_Type, Org_Type_Specific) %>% 
  unique() %>% mutate(Match = "match") %>%
  right_join(mfs.with.NA.in.orgs, by = "Match")%>%
  mutate(cluster_letters = "not observed", 
         std_area = NA) 

orgs.dat.long.med <- bind_rows(orgs.dat.long.med %>% ungroup %>% mutate(cluster_letters = as.character(cluster_letters)), matcher)

orgs.dat.long.med$cluster_letters <- factor(orgs.dat.long.med$cluster_letters, 
                                    levels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "not observed"))
orgs.dat.long.med$MassFeature_Column <- factor(orgs.dat.long.med$MassFeature_Column, 
                                               levels = rev(Org.MF.order$MassFeature_Column))


#Assign order of org samps-----
order.of.org.df <- Org.dat.long %>% select(SampID) %>% 
   left_join(meta.culture.dat %>% rename(SampID = CultureID), by = "SampID") %>%
   select(Org_Type, Org_Type_Specific, CultureID_short) %>%
   unique() %>% arrange(Org_Type, Org_Type_Specific)
order.of.org <- order.of.org.df$CultureID_short
order.of.org.df$CultureID_short <- factor(order.of.org.df$CultureID_short, 
                                          levels = (order.of.org))


orgs.dat.long.med$CultureID_short <- factor(orgs.dat.long.med$CultureID_short, 
                                               levels = (order.of.org))


#Plot up org tile plot-----
org.cluster.labs <- c("a", "b", "c", "d", "e", "f", "g", "")
names(org.cluster.labs) <- c("a", "b", "c", "d", "e", "f", "g", "not observed")

not.observe.label <-  data.frame(text.x = 18,text.y = 5,lab = "not observed", cluster_letters = factor("not observed",levels = names(org.cluster.labs)))

g.tile.org <- ggplot(data = orgs.dat.long.med %>% filter(CultureID_short != "Nmar"), aes(x = factor(CultureID_short), y = MassFeature_Column, fill = std_area, colour = "")) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  facet_grid(cluster_letters ~ ., scales="free_y", space="free_y", labeller = labeller(cluster_letters = org.cluster.labs))+
  geom_hline(data= lindat, aes(yintercept=x), col="black", size = 0.5) +
  geom_text(data = not.observe.label, inherit.aes = FALSE, aes(x = text.x, y =text.y, label = lab), fontface = "italic", size = 2)+
  scale_fill_gradient2(low="slateblue4", mid="cornsilk1", high="#00295D",
                       midpoint=0, limits = c(0.0, 1), na.value = "grey90", breaks = c(0.0, 0.5, 1.0))+
  scale_colour_manual(values=c("grey80")) +  
  guides(colour=guide_legend("not \nobserved", override.aes=list(fill="grey90")))+
  theme(axis.line.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
       # strip.text.placement = "outside",
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 8, face = "italic", angle=0),
        panel.spacing.y=unit(0, "lines"),
        legend.position = "top",
        legend.justification="right",
        legend.margin=margin(0,0,-15,0),
        legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        plot.margin = margin(.25, .25, 0, -.1, "cm"))+
  labs(x ="Organism", y = "Mass Feature", fill = "Standardized \n peak area") 
g.tile.org


#Make another tile plot to smash together with colored names-----
pal <- park_palette("Redwoods", 5)
#pal2 <- c("deepskyblue4", pal[5], pal[2:4], pal[1])
g.tile.org.key <- ggplot(data = order.of.org.df %>% filter(CultureID_short != "Nmar"), 
                         aes(x = factor(CultureID_short), y = 1, fill = Org_Type), colour = "black") +
  geom_tile() +
  geom_text(aes(label = CultureID_short), angle = -90, size = 1.7, fontface = "italic")+
  scale_fill_manual(values = pal)+
  xlab("Organisms (colored by broad classification)")+
  theme(axis.line.y = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = -4, r = 0, b = 0, l = 0),size = 7),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none",
        plot.margin = margin(-.1, 0, .5, -.1, "cm"))

g.tile.org.key


#Try to put the plots together!-----
g.clus.2 <- gclus + theme(plot.margin = margin(0, -.5, 1.6, -.2, "cm")) #trbl change bottom to match better!

g.tile.org.combo <- plot_grid(g.tile.org, g.tile.org.key, 
                              align = "v", axis = "lr",
                              rel_heights = c(1,0.11),
                              ncol = 1, scale = 1)
g.left.two <- plot_grid(g.tile.KOK, g.clus.2, 
                        align = "h", axis = "bt",
                        rel_widths = c(2, 0.7),
                        ncol = 2, scale = 1)
g.all <- plot_grid(g.left.two, g.tile.org.combo,
                  # align = "h", axis = "t",
                   rel_widths = c(3, 3),
                   rel_heights =  c(1, 2),
                   ncol = 2, scale = 1)
g.all


save_plot("Figures/Manuscript_figures/Clusters_and_Alluvial2.pdf", g.all, base_height = 16, base_width = 16, units="cm")











library(tidyverse)
library(scales)
library(cowplot)
library(here)
library(ggalluvial)



#Name your inputs
KOK.dat.wide.std.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
org.dat.wide.std.file <- "Intermediates/organs_wide_stand_withclusters.csv"
org.dat.dat <- "Intermediates/culture_combined_long.csv"
meta.dat.culture.filename <- "MetaData/CultureMetaData.csv"
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"

#######

#Read in your files
KOK.dat.long <- read_csv(KOK.dat.wide.std.file)
orgs.dat.long <- read_csv(org.dat.wide.std.file)
meta.culture.dat <- read_csv(meta.dat.culture.filename)

#get order of KOK samps
orderofSamps <- KOK.dat.long %>% select(SampID, latitude) %>% 
  unique() %>% 
  arrange(latitude) 
orderofSamps <- orderofSamps$SampID
KOK.dat.long$SampID <- factor(KOK.dat.long$SampID, levels = orderofSamps)


#get order of org samps
order.of.orgs <- orgs.dat.long %>% select(SampID) %>% 
  left_join(meta.culture.dat %>% rename(SampID = CultureID)) %>%
  select(SampID, Org_Type, Org_Type_Specific) %>%
  unique() %>% arrange(Org_Type, Org_Type_Specific)
order.of.orgs <- order.of.orgs$SampID

#get long org dat
org.dat.long.biovol <- read_csv(org.dat.dat) %>%
  rename(SampID = ID_rep) %>%
  select(MassFeature_Column, SampID, LogBioArea)

orgs.dat.long <- orgs.dat.long %>%
  left_join(org.dat.long.biovol) %>%
  mutate(LogBioArea = ifelse(LogBioArea == 0, NA, LogBioArea))
orgs.dat.long$SampID <- factor(orgs.dat.long$SampID, levels = order.of.orgs)

clusters.combined.KOK.org <- KOK.dat.long %>% select(MassFeature_Column, cluster) %>%
  rename(cluster.KOK = cluster) %>% unique() %>%
  left_join(orgs.dat.long %>% select(MassFeature_Column, cluster) %>%
              rename(cluster.org = cluster) %>% unique()) %>%
  group_by(cluster.KOK, cluster.org) %>%
  summarise(Freq = n()) %>% ungroup() %>%
  mutate(cluster.org = ifelse(is.na(cluster.org), "not observed", cluster.org))

cluster.summary.KOK <- KOK.dat.long %>% select(MassFeature_Column, cluster) %>%
  rename(cluster.KOK = cluster) %>% unique() %>%
  group_by(cluster.KOK) %>%
  summarise(total.KOK = n()) %>% ungroup() 

cluster.summary.org <- orgs.dat.long %>% select(MassFeature_Column, cluster)%>%
  mutate(cluster.org = as.character(cluster)) %>% unique()%>%
  group_by(cluster.org) %>%
  summarise(total.org = n()) %>% ungroup() 

clusters.combined.KOK.org <- clusters.combined.KOK.org %>%
  left_join(cluster.summary.KOK) %>%
  left_join(cluster.summary.org) %>%
  mutate(cluster.org = ifelse(is.na(cluster.org), "not observed", cluster.org))%>%
  mutate(total.org = ifelse(is.na(total.org), 55 , total.org)) %>%
  mutate(expected_Freq = total.org*total.KOK/sum(Freq))%>%
  mutate(higher_than_expected = ifelse(Freq>expected_Freq*1.2, 1, 0)) %>%
  mutate(higher_than_expected = ifelse(is.na(expected_Freq), -1, higher_than_expected) %>% as.character()) %>%
  mutate(cluster.KOK = paste0("KOK", cluster.KOK))

order.of.MF <-KOK.dat.long %>% select(MassFeature_Column, cluster) %>%
  rename(cluster.KOK = cluster) %>% unique() %>%
  left_join(orgs.dat.long %>% select(MassFeature_Column, cluster) %>%
              rename(cluster.org = cluster) ) %>% unique()

order.of.MF$cluster.org <- factor(order.of.MF$cluster.org, 
                                                levels = c("7","8","5","2","9","1","4", "6","3","10", "not observed"))

order.of.MF.KOK <- order.of.MF %>% arrange(cluster.KOK, cluster.org)
order.of.MF.org <- order.of.MF %>% arrange(cluster.org, cluster.KOK)
KOK.dat.long$MassFeature_Column <- factor(KOK.dat.long$MassFeature_Column, 
                       levels = rev(order.of.MF.org$MassFeature_Column))

#Make a tile plot of the clusters by KOK-------
g2 <- ggplot(dat = KOK.dat.long %>%
               mutate(std_area = ifelse(std_area > 3, 3, std_area)) %>%
               mutate(std_area = ifelse(std_area < -3, -3, std_area)), aes(x = SampID, y = MassFeature_Column, fill = std_area)) +
  geom_tile() +
  geom_tile(colour="white",size=0.1)+
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  scale_fill_gradient2(low="slateblue4", mid="white", high="firebrick2",
                       midpoint=0, limits = c(-3, 3), na.value = "grey")  +
  theme(axis.text.x = element_text(angle=-60, hjust=0),
        axis.ticks = element_blank(), axis.text.y = element_blank()) +
  theme(legend.position="bottom") +
  labs(x ="South <-> North", y = "Mass Feature")

#Make a tile plot of the clusters by organism (median reps) -----
orgs.dat.long.med <- orgs.dat.long %>% 
  left_join(meta.culture.dat %>% rename(SampID = CultureID)) %>%
  group_by(CultureID_short, MassFeature_Column, Org_Type, Org_Type_Specific, cluster) %>%
  mutate(LogBioArea = ifelse(is.na(LogBioArea), 0, LogBioArea)) %>%
  summarize(LogBioArea = median(LogBioArea)) %>%
  mutate(LogBioArea = ifelse(LogBioArea == 0, NA, LogBioArea))

order.of.orgs.med <- orgs.dat.long.med %>% ungroup() %>%
  select(CultureID_short, Org_Type, Org_Type_Specific) %>%
  unique() %>% arrange(Org_Type, Org_Type_Specific)
order.of.orgs.med <- order.of.orgs.med$CultureID_short
orgs.dat.long.med$CultureID_short <- factor(orgs.dat.long.med$CultureID_short, levels = order.of.orgs.med)

#Add in MFs that we didn't see in oganisms
mfs.with.NA.in.orgs <- KOK.dat.long %>% select(MassFeature_Column, cluster) %>%
  rename(cluster.KOK = cluster) %>% unique() %>%
  left_join(orgs.dat.long %>% select(MassFeature_Column, cluster) %>%
              rename(cluster.org = cluster)) %>%
  filter(is.na(cluster.org)) %>% unique() %>%
  mutate(Match = "match")

matcher <- orgs.dat.long.med %>% ungroup() %>%
  select(CultureID_short, Org_Type, Org_Type_Specific) %>% 
  unique() %>% mutate(Match = "match") %>%
  right_join(mfs.with.NA.in.orgs)%>%
  mutate(cluster = "not observed", 
         LogBioArea = NA) 

orgs.dat.long.med <- bind_rows(orgs.dat.long.med %>% mutate(cluster = as.character(cluster)), matcher)
orgs.dat.long.med$cluster <- factor(orgs.dat.long.med$cluster, 
                                                levels = c("7","8","5","2","9","1","4", "6","3","10", "not observed"))
orgs.dat.long.med$MassFeature_Column <- factor(orgs.dat.long.med$MassFeature_Column, 
                                               levels = rev(order.of.MF.org$MassFeature_Column))
orgs.dat.long.med$CultureID_short <- factor(orgs.dat.long.med$CultureID_short, levels = order.of.orgs.med)

lindat <- tibble(f=as.factor("1"), x=c(-1,1)*Inf)

g1.med <- ggplot(dat = orgs.dat.long.med, aes(x = CultureID_short, y = MassFeature_Column, fill = LogBioArea)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  geom_hline(data= lindat, aes(yintercept=x), col="black", size = 0.5) +
  scale_fill_gradient2(low="purple1", high="purple4", na.value = "gray95")+
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 7),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        panel.spacing.y=unit(0, "lines"),
        legend.position="bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.margin = margin(2, 0, 2, -.1, "cm"))+
  labs(x ="Organisms", y = "Mass Feature", fill = "Log(biovolume \n normalized \n peak area)")

#Make a tile plot of the KOK clusters (median reps) -----
KOK.dat.long.med <- KOK.dat.long %>% 
  mutate(latitude = as.factor(round(latitude, digits =1))) %>%
  group_by(latitude, MassFeature_Column, cluster) %>%
  summarize(std_area = median(std_area))

KOK.dat.long.med$MassFeature_Column <- factor(KOK.dat.long.med$MassFeature_Column, 
                                          levels = rev(order.of.MF.org$MassFeature_Column))

g2.med <- ggplot(dat = KOK.dat.long.med %>% 
               mutate(std_area = ifelse(std_area > 2.5, 2.5, std_area)) %>%
               mutate(std_area = ifelse(std_area < -2.5, -2.5, std_area)), aes(x = latitude, y = MassFeature_Column, fill = std_area)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  geom_hline(data= lindat, aes(yintercept=x), col="black", size = 0.5) +
  scale_fill_gradient2(low="slateblue4", mid="white", high="firebrick2",
                       midpoint=0, limits = c(-3, 3), na.value = "grey")  +
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 7),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_blank(),
        panel.spacing.y=unit(0, "lines"),
        legend.position="bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.margin = margin(2, 0, 2, 0, "cm"))+
  labs(x ="Latitude", y = "Mass Feature", fill = "Standardized \n peak area")


#Make an alluvial plot of organisms and KOK clusters----

clusters.combined.KOK.org$cluster.org <- factor(clusters.combined.KOK.org$cluster.org, 
                                                   levels = c("7","8","5","2","9","1","4", "6","3","10", "not observed"))


gclus<- ggplot(data = clusters.combined.KOK.org,
           aes(y = Freq, axis1 = cluster.KOK, axis2 = cluster.org)) +
  geom_alluvium(aes(alpha = higher_than_expected), fill = "black", width = 1/12,  lode.guidance = "forward") +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  scale_x_discrete(limits = c("KOK1606_clusters", "Organism_clusters"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
    scale_alpha_manual(values = c("0.2",  "1"))+
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
    legend.position = "none",
    plot.margin = margin(2, 0, 2, 0, "cm"))


#Combine them
g.all <- plot_grid(g2.med, gclus, g1.med,
                   align = "h", axis = "bt",
                   rel_widths = c(2, 1, 2),
                   ncol = 3, scale = 1)

save_plot("Figures/Preliminary/Clusters_and_Alluvial.pdf", g.all, base_height = 7.5, base_aspect_ratio = 0.8)

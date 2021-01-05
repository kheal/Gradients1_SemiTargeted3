library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(RCurl)

#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         BestMatch = Compound.Name_figure) %>%
  select(BestMatch, Identification) %>% unique()

#Load data----
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
MF.dat.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
Quan.dat.file <- "Intermediates/Quantified_MFSummary_Enviro.csv"

#Get MGL cluster info----
MGL.dat <- read_csv(MGL.dat.file) 

MGL.dat <- MGL.dat %>%
  mutate(MGL_clusters = paste0("MGL_", cluster_letters)) %>% 
  select(MassFeature_Column, MGL_clusters) %>%
  unique() 

#Get KM cluster info----
KM.dat <- read_csv(KM.dat.file)

KM.dat <- KM.dat %>%
  mutate(KM_clusters = paste0("KM_", cluster_letters)) %>% 
  select(MassFeature_Column, KM_clusters) %>%
  unique()

#Get KOK cluster info----
KOK.dat <- read_csv(KOK.dat.file)

KOK.dat <- KOK.dat %>%
  mutate(KOK_clusters = paste0("KOK_", cluster_letters)) %>% 
  select(MassFeature_Column, KOK_clusters) %>%
  unique()

#Load MF information----
MF.dat <- read_csv(MF.dat.file) %>%
  filter(Confidence == 1) %>% 
  select(MassFeature_Column, Identification) %>%  unique() %>%
  left_join(KOK.dat) %>% left_join(KM.dat) %>% left_join(MGL.dat) %>%
  left_join(stds.dat, by = "Identification") %>% 
  select(-Identification) %>%
  rename(Identification = BestMatch)

#Get QuanDat info on here so we can pull out the high abundance ones that fit each pattern----
dat.quan <- read_csv(Quan.dat.file) %>%
  select(Identification, nmolCmed) %>%
  left_join(stds.dat, by = "Identification") %>% 
  select(-Identification) %>%
  rename(Identification = BestMatch)
MF.dat.quan <- MF.dat %>% left_join(dat.quan) %>%
  arrange(desc(nmolCmed))
MF.dat.quan$Identification <- factor(MF.dat.quan$Identification, levels = MF.dat.quan$Identification)


#First for the KOK data----
MF.dat.quan.KOK <- MF.dat.quan %>%
  select(Identification, KOK_clusters) %>%
  arrange(KOK_clusters, Identification) %>%
  mutate(Identification = factor(Identification, levels = Identification)) %>%
  group_by(KOK_clusters) %>%
  mutate(row.plot = order(Identification)) %>% ungroup()

MF.dat.quan.KOK <- rbind(MF.dat.quan.KOK) %>%
  mutate(KOK_clusters = str_remove(KOK_clusters, "KOK_")) %>%
  mutate(column.plot = ifelse(row_number() > 33, 1, 1)) %>%
  mutate(column.plot.real = ifelse(KOK_clusters == "a", column.plot, 1)) %>%
  mutate(row.plot = ifelse(column.plot.real == 1, row_number()-33, row.plot)) 

g.cmps.KOK <- ggplot() +
  geom_text(data = MF.dat.quan.KOK, aes(y = fct_rev(Identification), x = 1, label = Identification), size = 1.8) +
  facet_wrap(KOK_clusters ~ ., ncol = 7, scales = "free") +
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
g.cmps.KOK

#Next for the MGL data----
MF.dat.quan.MGL <- MF.dat.quan %>%
  select(Identification, MGL_clusters)
MF.dat.quan.MGL.extra1 <- c(NA, NA)
MF.dat.quan.MGL.extra2 <- c(NA, NA)
MF.dat.quan.MGL <- rbind(MF.dat.quan.MGL, MF.dat.quan.MGL.extra1) %>%
  rbind(MF.dat.quan.MGL.extra2) %>%
  mutate(MGL_clusters = str_remove(MGL_clusters, "MGL_")) %>%
  filter(!is.na(MGL_clusters))

g.cmps.MGL <- ggplot() +
  geom_text(data = MF.dat.quan.MGL, aes(y = fct_rev(Identification), x = 1, label = Identification), size = 1.8) +
  facet_wrap(MGL_clusters ~ ., ncol = 7, scales = "free") +
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
g.cmps.MGL


#Finally for the KM data----
MF.dat.quan.KM <- MF.dat.quan %>%
  select(Identification, KM_clusters) %>%
  arrange(KM_clusters, Identification) %>%
  mutate(Identification = factor(Identification, levels = Identification)) %>%
  group_by(KM_clusters) %>%
  mutate(row.plot = order(Identification)) %>% ungroup() %>%
  filter(!is.na(KM_clusters))

MF.dat.quan.KM <- rbind(MF.dat.quan.KM) %>%
  mutate(KM_clusters = str_remove(KM_clusters, "KM_")) %>%
  mutate(column.plot = ifelse(row_number() > 51, 1.5, 1)) %>%
  mutate(column.plot.real = ifelse(KM_clusters == "c", column.plot, 1)) %>%
  mutate(row.plot = ifelse(column.plot.real == 1.5, row_number()-51, row.plot)) 

g.cmps.KM <- ggplot() +
  geom_text(data = MF.dat.quan.KM, aes(y = fct_rev(Identification), x = 1, label = Identification), size = 1.8) +
  facet_wrap(KM_clusters ~ ., ncol = 7, scales = "free") +
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
g.cmps.KM


rm(list=setdiff(ls(), c("g.mode.KOK", "g.mode.KM", "g.mode.MGL", "g.MGL.ctd", "g.KM.ctd", "g.tile.KM", "g.tile.MGL", "g.tile.KOK", "g.KOK.ctd", "g.cmps.KOK", "g.cmps.MGL", "g.cmps.KM")))


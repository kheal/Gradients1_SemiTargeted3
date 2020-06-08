#Helpfulwebsite: https://briatte.github.io/ggnet/

library(tidyverse)
library(cowplot)
library(here)
options(readr.num_columns = 0)

#File names
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Org.dat.file <- "Intermediates/organs_wide_stand_withclusters.csv"

MGL.KM.boot.file <- "Intermediates/BootstrapResults/KMvsMGL.csv"
MGL.KOK.boot.file <- "Intermediates/BootstrapResults/KOKvsMGL.csv"
MGL.Org.boot.file <- 
KM.KOK.boot.file <- "Intermediates/BootstrapResults/KOKvsKM.csv"
KM.Org.boot.file <- 
KM.KOK.boot.file <- "Intermediates/BootstrapResults/KOKvsKM.csv"



#Load files
MGL.dat <- read_csv(MGL.dat.file) %>% mutate(cluster.letters.MGL = paste0(cluster_letters, "_MGL")) %>% 
  select(MassFeature_Column, cluster.letters.MGL) %>% unique()
KM.dat <- read_csv(KM.dat.file) %>% mutate(cluster.letters.KM = paste0(cluster_letters, "_KM"))%>% 
  select(MassFeature_Column, cluster.letters.KM) %>% unique()
KOK.dat <- read_csv(KOK.dat.file) %>% mutate(cluster.letters.KOK = paste0(cluster_letters, "_KOK"))%>% 
  select(MassFeature_Column, cluster.letters.KOK) %>% unique()
Org.dat <- read_csv(Org.dat.file) %>% mutate(cluster.letters.Orgs = paste0(cluster_letters, "_Orgs"))%>% 
  select(MassFeature_Column, cluster.letters.Orgs) %>% unique()

All.clusters <- KOK.dat %>% left_join(MGL.dat, by = "MassFeature_Column") %>% 
  left_join(KM.dat, by = "MassFeature_Column") %>% 
  left_join(Org.dat, by = "MassFeature_Column") %>%
  mutate(cluster.letters.MGL = ifelse(is.na(cluster.letters.MGL), "NA_MGL", cluster.letters.MGL)) %>%
  mutate(cluster.letters.Orgs = ifelse(is.na(cluster.letters.Orgs), "NA_Orgs", cluster.letters.Orgs)) %>%
  mutate(cluster.letters.KM = ifelse(is.na(cluster.letters.KM), "NA_KM", cluster.letters.KM))

#This will be our verticies information 
Cluster.info <- All.clusters %>%
  gather(key = sample.set, value = cluster, -MassFeature_Column) %>%
  group_by(sample.set, cluster) %>%
  summarise(count = n())

#Now we have to go and get the connections from the bootstrapping information
MGL.KM.boot <- read_csv(MGL.KM.boot.file)

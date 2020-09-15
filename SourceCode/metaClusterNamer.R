#This code is to assign each mass feature into the 3 Meta-clusters (if applicable)
library(tidyverse)

#Name in your dat files
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Org.dat.file <- "Intermediates/organs_wide_stand_withclusters.csv"
MF.dat.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"

#Read in your dat files
MGL.dat <- read_csv(MGL.dat.file) %>% select(MassFeature_Column, cluster_letters) %>%
  rename(MGL_clusters = cluster_letters) %>% unique() 
KM.dat <- read_csv(KM.dat.file) %>% select(MassFeature_Column, cluster_letters) %>%
  rename(KM_clusters = cluster_letters) %>% unique()
KOK.dat <- read_csv(KOK.dat.file) %>% select(MassFeature_Column, cluster_letters) %>%
  rename(KOK_clusters = cluster_letters) %>% unique()
Org.dat <- read_csv(Org.dat.file) %>% select(MassFeature_Column, cluster_letters) %>%
  rename(Org_clusters = cluster_letters) %>% unique()
MF.dat <- read_csv(MF.dat.file)

#Smash them together
combo.dat <- KOK.dat %>%
  left_join(MGL.dat, by = "MassFeature_Column") %>%
  left_join(KM.dat, by = "MassFeature_Column") %>%
  left_join(Org.dat, by = "MassFeature_Column")

#Make some rules
combo.dat.rare <- combo.dat %>%
  mutate(Rare_org = ifelse(is.na(Org_clusters), 1, 0),
         Rare_KOK = ifelse(KOK_clusters == "f" | KOK_clusters == "g", 1, 0),
         Rare_MGL = ifelse(is.na(MGL_clusters) | MGL_clusters == "f", 1, 0),
         Rare_KM = ifelse(is.na(KM_clusters) | KM_clusters == "d", 1, 0)) %>%
  mutate(Rare_field = Rare_KOK+Rare_MGL+Rare_KM) %>%
  mutate(MetaClusterAssignment = ifelse(Rare_field > 1 & Rare_org == 1, "Rare", NA)) %>%
  select(MassFeature_Column, MetaClusterAssignment) %>%
  filter(!is.na(MetaClusterAssignment))

combo.dat.core <- combo.dat %>%
  mutate(Core_org = ifelse(Org_clusters =="a", 1, 0),
         Core_KOK = ifelse(KOK_clusters == "a", 1, 0),
         Core_MGL = ifelse(MGL_clusters == "b", 1, 0),
         Core_KM = ifelse(KM_clusters == "a", 1, 0)) %>%
  mutate(Core_field = Core_KOK+Core_MGL+Core_KM) %>%
  mutate(MetaClusterAssignment = ifelse(Core_field > 1 & Core_org == 1, "Core", NA))%>%
  select(MassFeature_Column, MetaClusterAssignment) %>%
  filter(!is.na(MetaClusterAssignment))


combo.dat.dino <- combo.dat %>%
  mutate(Dino_org = ifelse(Org_clusters =="e", 1, 0),
         Dino_KOK = ifelse(KOK_clusters == "b" | KOK_clusters == "c", 1, 0),
         Dino_MGL = ifelse(MGL_clusters == "c", 1, 0),
         Dino_KM = ifelse(KM_clusters == "c", 1, 0)) %>%
  mutate(Dino_field = Dino_KOK+Dino_MGL+Dino_KM) %>%
  mutate(MetaClusterAssignment = ifelse(Dino_field > 1 & Dino_org == 1, "Dino", NA))%>%
  select(MassFeature_Column, MetaClusterAssignment) %>%
  filter(!is.na(MetaClusterAssignment))

#Clean up
combo.dat.metaclus <- combo.dat %>%
  left_join(bind_rows(combo.dat.rare, bind_rows(combo.dat.core, combo.dat.dino)))

#Write it out
write_csv(combo.dat.metaclus, "Intermediates/metacluster_assignments.csv")

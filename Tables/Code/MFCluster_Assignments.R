#TO DO: make table that is the known compounds and which assignments they live in.
library(tidyverse)

#Read in your dat files
dat.filename <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
MGL.clusterfilename  <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.clusterfilename <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.clusterfilename <- "Intermediates/KOK_wide_stand_withclusters.csv"
Org.clusterfilename <- "Intermediates/organs_wide_stand_withclusters.csv"

#Check out data
dat <- read_csv(dat.filename) %>%
  select(MassFeature_Column, Identification:z)
MGL.dat <- read_csv(MGL.clusterfilename) %>%
  select(MassFeature_Column, cluster_letters) %>% unique() %>%
  rename(`NPTZ depth profile cluster` = cluster_letters)
KM.dat <- read_csv(KM.clusterfilename) %>%
  select(MassFeature_Column, cluster_letters) %>% unique() %>%
  rename(`NPSG depth profile cluster` = cluster_letters)
KOK.dat <- read_csv(KOK.clusterfilename) %>%
  select(MassFeature_Column, cluster_letters) %>% unique()%>%
  rename(`Transect cluster` = cluster_letters)
Org.dat <- read_csv(Org.clusterfilename) %>%
  select(MassFeature_Column, cluster_letters) %>% unique()%>%
  rename(`Culture cluster` = cluster_letters)

dat.combo <- dat %>%
  left_join(MGL.dat, by = "MassFeature_Column") %>%
  left_join(KM.dat, by = "MassFeature_Column") %>%
  left_join(KOK.dat, by = "MassFeature_Column") %>%
  left_join(Org.dat, by = "MassFeature_Column") %>%
  arrange(mz) %>%
  arrange(Column) %>%
  arrange(z)
  

write_csv(dat.combo, "Tables/Manuscript_tables/SuppTables/MFCluster_Assignments.csv")
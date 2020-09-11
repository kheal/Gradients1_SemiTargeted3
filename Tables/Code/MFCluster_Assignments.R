library(tidyverse)

#Read in your dat files
dat.filename <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
dat.ms2.filename <- "Intermediates/MF_info_withMS2s.csv"
MGL.clusterfilename  <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.clusterfilename <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.clusterfilename <- "Intermediates/KOK_wide_stand_withclusters.csv"
Org.clusterfilename <- "Intermediates/organs_wide_stand_withclusters.csv"
Metacluster.filename <- "Intermediates/metacluster_assignments.csv"

#Check out MS2 data
dat.ms2 <- read_csv(dat.ms2.filename)

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
metacluster.dat <- read_csv(Metacluster.filename)

dat.combo <- dat %>%
  left_join(MGL.dat, by = "MassFeature_Column") %>%
  left_join(KM.dat, by = "MassFeature_Column") %>%
  left_join(KOK.dat, by = "MassFeature_Column") %>%
  left_join(Org.dat, by = "MassFeature_Column") %>%
  left_join(metacluster.dat, by = "MassFeature_Column")%>%
  arrange(mz) %>%
  arrange(Column) %>%
  arrange(z)

dat.combo.MS2 <- dat.combo %>%
  left_join(dat.ms2, by = c("MassFeature_Column", "Identification", "Confidence", "mz", "rt", "Column", "z"))

write_csv(dat.combo.MS2, "Tables/Manuscript_tables/SuppTables/MFCluster_Assignments.csv")

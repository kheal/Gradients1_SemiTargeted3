library(here)
library(tidyverse)

#Name data----
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Org.dat.file <- "Intermediates/organs_wide_stand_withclusters.csv"
ID.dat.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"

#Load data----
KM.dat <- read_csv(KM.dat.file)
ID.dat <- read_csv(ID.dat.file)
MGL.dat <- read_csv(MGL.dat.file)
KOK.dat <- read_csv(KOK.dat.file)

#Mudge data for KM ----
ID.dat.sub <- ID.dat %>% select(MassFeature_Column, mz, rt)

KM.dat.wide <- KM.dat %>%
  select(SampID, std_area, MassFeature_Column) %>%
  spread(data = ., key = SampID, value = std_area) %>%
  left_join(ID.dat.sub) 

KM.dat.wide.HILICNeg <- KM.dat.wide %>%
  filter(str_detect(MassFeature_Column, "Neg")) %>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())

KM.dat.wide.HILICPos <- KM.dat.wide %>%
  filter(str_detect(MassFeature_Column, "Pos"))%>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())

KM.dat.wide.RP <- KM.dat.wide %>%
  filter(str_detect(MassFeature_Column, "RP"))%>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())


#Mudge data for KOK ----
ID.dat.sub <- ID.dat %>% select(MassFeature_Column, mz, rt)

KOK.dat.wide <- KOK.dat %>%
  select(SampID, std_area, MassFeature_Column) %>%
  spread(data = ., key = SampID, value = std_area) %>%
  left_join(ID.dat.sub) 

KOK.dat.wide.HILICNeg <- KOK.dat.wide %>%
  filter(str_detect(MassFeature_Column, "Neg")) %>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())

KOK.dat.wide.HILICPos <- KOK.dat.wide %>%
  filter(str_detect(MassFeature_Column, "Pos"))%>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())

KOK.dat.wide.RP <- KOK.dat.wide %>%
  filter(str_detect(MassFeature_Column, "RP"))%>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())


#Mudge data for MGL ----
ID.dat.sub <- ID.dat %>% select(MassFeature_Column, mz, rt)

MGL.dat.wide <- MGL.dat %>%
  select(SampID, std_area, MassFeature_Column) %>%
  spread(data = ., key = SampID, value = std_area) %>%
  left_join(ID.dat.sub) 

MGL.dat.wide.HILICNeg <- MGL.dat.wide %>%
  filter(str_detect(MassFeature_Column, "Neg")) %>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())

MGL.dat.wide.HILICPos <- MGL.dat.wide %>%
  filter(str_detect(MassFeature_Column, "Pos"))%>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())

MGL.dat.wide.RP <- MGL.dat.wide %>%
  filter(str_detect(MassFeature_Column, "RP"))%>%
  mutate(mz_rt = paste0(mz, "_", rt)) %>%
  select(-MassFeature_Column, -mz,-rt) %>%
  select(mz_rt, everything())


#Save out data for KM ----
write_csv(KM.dat.wide.HILICNeg, "Intermediates/Metabolomics_Workbench_processedData/KM1513_HILICNeg")
write_csv(KM.dat.wide.HILICPos, "Intermediates/Metabolomics_Workbench_processedData/KM1513_HILICPos")
write_csv(KM.dat.wide.RP, "Intermediates/Metabolomics_Workbench_processedData/KM1513_HILICRP")

#Save out data for KOK ----
write_csv(KOK.dat.wide.HILICNeg, "Intermediates/Metabolomics_Workbench_processedData/KOK1606_HILICNeg")
write_csv(KOK.dat.wide.HILICPos, "Intermediates/Metabolomics_Workbench_processedData/KOK1606_HILICPos")
write_csv(KOK.dat.wide.RP, "Intermediates/Metabolomics_Workbench_processedData/KOK1606_HILICRP")

#Save out data for MGL ----
write_csv(MGL.dat.wide.HILICNeg, "Intermediates/Metabolomics_Workbench_processedData/MGL1704_HILICNeg")
write_csv(MGL.dat.wide.HILICPos, "Intermediates/Metabolomics_Workbench_processedData/MGL1704_HILICPos")
write_csv(MGL.dat.wide.RP, "Intermediates/Metabolomics_Workbench_processedData/MGL1704_HILICRP")






#Get list of problem mass features to go and double check the integrations on-----
ID.problem.masses <- ID.dat %>%
  mutate(mz.round = round(mz, digits = 2)) %>%
  group_by(mz.round, Column, z) %>%
  summarise(count = n()) %>%
  filter(count > 1)

Problem.MFs <- ID.dat %>%
  mutate(mz.round = round(mz, digits = 2)) %>%
  filter(mz.round %in% ID.problem.masses$mz.round) %>%
  mutate(RT = rt/60) %>%
  select(MassFeature_Column, Identification, Column, z, mz, RT, rt) %>%
  arrange(Column, z, mz)

write_csv(Problem.MFs, "Intermediates/problemMFs_tocheck.csv")

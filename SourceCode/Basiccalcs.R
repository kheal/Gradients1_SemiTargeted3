#This code performs basic calculations for the results section of the manuscript programatically
library(tidyverse)
library(fuzzyjoin)

#Load data----
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Org.dat.file <- "Intermediates/organs_wide_stand_withclusters.csv"
Quan.dat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
Quan.dat.file <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
Meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
meta.dat.culture.file <- "MetaData/CultureMetaData.csv"
PC.dat.file <- "MetaData/PCPN/KOK1606_PCPN_UW_Preliminary_OSU_KRH.csv"

#Number of "quality mass features"----
kok.dat <- read_csv(KOK.dat.file) 
MF.count <- kok.dat %>% select(MassFeature_Column) %>%
  unique()
print(paste0("Total number quality mass features = ", length(MF.count$MassFeature_Column)))


#Number of compounds quantified----
quan.dat <- read_csv(Quan.dat.file) 
quan.compounds.count <- quan.dat %>% 
  select(Identification, nmolinEnviroave) %>% 
  mutate(quantified = ifelse(is.na(nmolinEnviroave), "no", "yes")) %>%
  filter(quantified == "yes") %>%
  select(-nmolinEnviroave) %>%
  unique()
print(paste0("Total number compounds ID = ", length(quan.compounds.count$Identification)))

#Range of metabolites quantified from each sample set----
quan.KOK <- quan.dat %>% 
  filter(!str_detect(SampID, "S7")) %>%
  filter(!str_detect(SampID, "KM1513")) %>%
  select(SampID, totalCmeasured_nM, totalNmeasured_nM) %>% unique()
print(paste0("Range of quantified metabolites from KOK cruise = ", round(min(quan.KOK$totalCmeasured_nM), digits = 1), " nM to ", round(max(quan.KOK$totalCmeasured_nM), digits = 1), " nM"))

quan.MGL <- quan.dat %>% 
  filter(str_detect(SampID, "S7")) %>%
  select(SampID, totalCmeasured_nM, totalNmeasured_nM) %>% unique()
print(paste0("Range of quantified metabolites from MGL samples = ", round(min(quan.MGL$totalCmeasured_nM), digits = 1), " to ", round(max(quan.MGL$totalCmeasured_nM), digits = 1), " nM"))

quan.KM <- quan.dat %>% 
  filter(str_detect(SampID, "KM1513")) %>%
  select(SampID, totalCmeasured_nM, totalNmeasured_nM) %>% unique()
print(paste0("Range of quantified metabolites from KM samples = ", round(min(quan.KM$totalCmeasured_nM), digits = 1), " nM to ", round(max(quan.KM$totalCmeasured_nM), digits = 1), " nM"))

#Percent of particulate carbon across the KOK gradient----
meta.dat <- read_csv(Meta.dat.file) %>%
  select(SampID, Station_1, latitude)

quan.dat.KOK.summ <- quan.KOK %>%
  left_join(meta.dat, by = "SampID") %>%
  group_by(Station_1) %>%
  summarise(nmolave.metabC = mean(totalCmeasured_nM),
            nmolstd.metabC = sd(totalCmeasured_nM), 
            nmolave.metabN = mean(totalNmeasured_nM),
            nmolstd.metabN = sd(totalNmeasured_nM), 
            latitude = mean(latitude)) %>%
  mutate(nmolstd.metabC = ifelse(is.na(nmolstd.metabC), 0, nmolstd.metabC)) %>%
  mutate(nmolstd.metabN = ifelse(is.na(nmolstd.metabN), 0, nmolstd.metabN))

PC.dat <- read_csv(PC.dat.file)

PC.dat.summary <- PC.dat %>%
  mutate(latitude.round = round(Latitude, digits = 0)) %>%
  group_by(latitude.round) %>%
  summarise(nmolave.PC = mean(PC)*1000,
            nmolstd.PC = sd(PC)*1000, 
            nmolave.PN = mean(PN)*1000,
            nmolstd.PN = sd(PN)*1000, 
            latitude = mean(Latitude))

mesh.KOK.PC <- quan.dat.KOK.summ %>%
  difference_left_join(PC.dat.summary, by = "latitude", max_dist = 0.5,
                       distance_col = "lat.difference") %>%
  filter(Station_1 != 7 | lat.difference < 0.47) %>%
  filter(Station_1 != 8 | lat.difference < 0.46)%>%
  mutate(fraction.PC = nmolave.metabC/nmolave.PC) %>%
  mutate(fraction.PC.sd = fraction.PC*(nmolstd.metabC/nmolave.PC+nmolstd.PC/nmolave.PC)) %>%
  mutate(percent.PC = fraction.PC*100, 
         sd.percent.PC = fraction.PC.sd*100) %>%
  mutate(fraction.PN = nmolave.metabN/nmolave.PN) %>%
  mutate(fraction.PN.sd = fraction.PN*(nmolstd.metabN/nmolave.PN+nmolstd.PN/nmolave.PN)) %>%
  mutate(percent.PN = fraction.PN*100, 
         sd.percent.PN = fraction.PN.sd*100) %>%
  select(Station_1, latitude.x, percent.PC, sd.percent.PC, percent.PN, sd.percent.PN)

g.pc <- ggplot() +
  geom_line(data = mesh.KOK.PC, aes(x = latitude.x, y =  percent.PC)) +
  geom_ribbon(data = mesh.KOK.PC, aes(x = latitude.x, ymin = percent.PC-sd.percent.PC, ymax = percent.PC+sd.percent.PC), alpha = 0.5)

print(g.pc)

#Percent of particulate nitrogen across the KOK gradient----
g.pn <- ggplot() +
  geom_line(data = mesh.KOK.PC, aes(x = latitude.x, y =  percent.PN)) +
  geom_ribbon(data = mesh.KOK.PC, aes(x = latitude.x, ymin = percent.PN-sd.percent.PN, ymax = percent.PN+sd.percent.PN), alpha = 0.5)

print(g.pn)

#Homarine calculations-----
meta.dat.culture <- read_csv(meta.dat.culture.file)
quan_dat_cultures <- read_csv(Quan.dat.file) %>%
  filter(Identification == "Homarine") %>%
 # select(Identification, ID_rep, Org_Name, Org_Type, C, Org_Type_Specific,nmolmetab_perC) %>%
  mutate(moleFractionHomarine = nmolmetab_perC/1000*100)

quan_dat_cultures2 <- quan_dat_cultures %>%
  rename(CultureID  = ID_rep) %>%
  left_join(meta.dat.culture, by = "CultureID") 





dat4cul <-  dat3cul %>%
  group_by(CultureID, Identification, Org_Type) %>%
  summarise(nmolmetab_perC_cul = mean(nmolmetab_perC, na.rm = T)) %>%
  group_by(Identification, Org_Type) %>%
  summarise(nmolmetab_perC_cul_med = median(nmolmetab_perC_cul, na.rm = T),
            nmolmetab_perC_cul_max = ifelse(is.na(nmolmetab_perC_cul_med), NA, max(nmolmetab_perC_cul, na.rm = T)), 
            nmolmetab_perC_cul_min = ifelse(is.na(nmolmetab_perC_cul_med), NA, min(nmolmetab_perC_cul, na.rm = T))) %>%
  filter(!is.na(nmolmetab_perC_cul_med)) 
dat5cul <-  dat3cul %>%
  group_by(CultureID, Identification) %>%
  summarise(nmolmetab_perC_cul = mean(nmolmetab_perC, na.rm = T)) %>%
  group_by(Identification) %>%
  summarise(nmolmetab_perC_cul_med = median(nmolmetab_perC_cul, na.rm = T),
            nmolmetab_perC_cul_max = ifelse(is.na(nmolmetab_perC_cul_med), NA, max(nmolmetab_perC_cul, na.rm = T)), 
            nmolmetab_perC_cul_min = ifelse(is.na(nmolmetab_perC_cul_med), NA, min(nmolmetab_perC_cul, na.rm = T))) %>%
  filter(!is.na(nmolmetab_perC_cul_med)) 

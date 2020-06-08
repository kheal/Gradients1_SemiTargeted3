#TO DO: Rerun and rewrite this so we're using the same RF and RFratios for the cultures and the environmental samples


library(tidyverse)
library(RCurl)
library(ggforce)
library(cowplot)

#Name your files
id.file <- "Intermediates/WideArea_withIDinfo.csv"
# id.manual.file <- "RawOutput/MFs_match_wManual.csv"
std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"
dat.hilicpos.file <- "Intermediates/Skyline_Standards/Standards_integrated_HILICPos.csv"
dat.hilicneg.file <-  "Intermediates/Skyline_Standards/Standards_integrated_HILICNeg.csv"
dat.cyano.file <- "Intermediates/Skyline_Standards/Standards_integrated_CyanoAq.csv"
dat.hilicpos.file.nostds <- "RawOutput/HILICPos_IntegrationsBigPeaksWTargeted2.csv"
dat.hilicneg.file.nostds <-  "RawOutput/HILICNeg_IntegrationsBigPeaksWTargeted.csv"
dat.cyano.file.nostds <- "RawOutput/CyanoAq_IntegrationsBigPeaksWTargeted.csv"
rf.ratio.matcher.file <- "Intermediates/RFs_relativeRF_matcher.csv"
samp.key.file <- "MetaData/Sample_key.csv"
dat.file <- "Intermediates/Longdata.csv"
is.names.file <- "MetaData/InternalStandardNames.csv"
meta.data.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"


#Get info on the standards that we can quantify -----
MF_info <- read.csv(id.file) %>%
  select(MassFeature_Column:z) %>%
  filter(Confidence == 1)

#Find which mix its in, what concentration its in that mix, add that to MF_info2----
IngallsStandards <- read.csv(text = getURL(std.url), header = T) 

## Use 1 = RF ratio for MTA since it's so late in the chromatogram----
IngallsStandards <- IngallsStandards %>%
  mutate(QE.RF.ratio = ifelse(Compound.Name == "Methylthioadenosine" | Compound.Name == "Trimethylammonium Propionate (TMAP)"
                              , 1, as.numeric(as.character(QE.RF.ratio))))
IngallsStandards2 <- IngallsStandards %>%
  mutate(Column = ifelse(Column == "RP", "RP", "HILIC"),
         Identification = Compound.Name) %>%
  select(Identification, Column, QE.RF.ratio, HILICMix, Conc..uM, Emperical.Formula, z)

MF_info2 <- MF_info %>%
  left_join(IngallsStandards2)%>%
  mutate(Conc..uM = ifelse(MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos", 4, Conc..uM))%>%
  filter(!is.na(Conc..uM)) %>%
  mutate(MassFeature_Column = MassFeature_Column %>% 
           str_replace_all("HILICNeg", "HILIC") %>%
           str_replace_all("HILICPos", "HILIC"))

#Get raw Skyline data of just the standard runs----
dat1 <- read.csv(dat.hilicpos.file) %>%
  mutate(Column = "HILIC", z = 1) %>% select(-Protein) 
dat2 <- read.csv(dat.hilicneg.file)  %>%
  mutate(Column = "HILIC", z = -1) %>% select(-Protein)
dat3 <- read.csv(dat.cyano.file)  %>%
  mutate(Column = "RP", z = 1) %>% select(-Protein)


#Meshing this values with the MF_info2
dat <- rbind(dat1, dat2) %>%
  bind_rows(dat3) %>%
  mutate(MF_Frac = Precursor.Ion.Name) %>%
  mutate(MassFeature_Column = paste(MF_Frac, Column, sep = "_X_"))%>%
  mutate(MassFeature_Column = MassFeature_Column %>% 
           str_replace_all("HILICNeg", "HILIC") %>%
           str_replace_all("HILICPos", "HILIC")) %>%
  right_join(MF_info2)

#Summarize RF for each compound by date ----  
dat.RF.bydate <- dat %>%
  mutate(date = str_extract(Replicate.Name, "^\\d{6}")) %>%
  mutate(HILICMix = as.character(HILICMix)) %>%
  mutate(HILICMix = ifelse(is.na(HILICMix), "RP", HILICMix)) %>%
  mutate(HILICMixMatch = str_extract(Replicate.Name, "Mix\\d")) %>%
  mutate(HILICMixMatch = ifelse(is.na(HILICMixMatch), HILICMix, HILICMixMatch)) %>%
  filter(!str_detect(Replicate.Name, "atrix")) %>%
  filter(HILICMix == as.character(HILICMixMatch)) %>%
  mutate(RF = as.numeric(Area)/as.numeric(Conc..uM)) %>%
  select(Column:RF) %>% 
  group_by(MassFeature_Column, date) %>%
  mutate(RF = mean(RF, na.rm = TRUE)) %>%
  unique()

write_csv(dat.RF.bydate, "Intermediates/RF_byDate_uncorrected.csv")

#Plot up these ----  SAVE THIS OUT
#mutate(RFmin = ifelse(MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos", 615089984, RFmin),
#       RFmax = ifelse(MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos", 615089984, RFmax))
Saveplot1 <- ggplot(data = dat.RF.bydate, aes(x = date, y = RF))+
  geom_bar(stat = "identity") +
  facet_wrap_paginate(vars(Identification), scales = "free",
                      ncol = 6, nrow = 6, page = 1)+
  theme(text = element_text(size=8), axis.text = element_text(size = 8))
Saveplot2 <- ggplot(data = dat.RF.bydate, aes(x = date, y = RF))+
  geom_bar(stat = "identity") +
  facet_wrap_paginate(vars(Identification), scales = "free",
                      ncol = 6, nrow = 6, page = 2)+
  theme(text = element_text(size=8), axis.text = element_text(size = 8))
Saveplot3 <- ggplot(data = dat.RF.bydate, aes(x = date, y = RF))+
  geom_bar(stat = "identity") +
  facet_wrap_paginate(vars(Identification), scales = "free",
                      ncol = 6, nrow = 6, page = 3)+
  theme(text = element_text(size=8), axis.text = element_text(size = 8))

dir.create(file.path("Intermediates", "RFplots"), showWarnings = FALSE)
save_plot(filename = "Intermediates/RFplots/RF_plot1.png",
       plot = Saveplot1,  base_width = 11, base_height = 8.5, units = "in", device = "png")
save_plot(filename = "Intermediates/RFplots/RF_plot2.png",
          plot = Saveplot2,  base_width = 11, base_height = 8.5, units = "in", device = "png")
save_plot(filename = "Intermediates/RFplots/RF_plot3.png",
          plot = Saveplot3,  base_width = 11, base_height = 8.5, units = "in", device = "png")

rm(Saveplot1, Saveplot2, Saveplot3)

#Get RF ratio for HILIC compounds - Use the data from 180702----
RFratio.HILIC <- dat %>%
  filter(!Column == "CyanoAq") %>%
  filter(str_detect(Replicate.Name, "180702|180709")) %>%
  filter(MassFeature_Column != "I132.102R9.3_HILIC_X_HILIC") %>%
  mutate(Type = paste(Standards = ifelse(str_detect(Replicate.Name, as.character(HILICMix)), "Standards", "Water"),     
                      Matrix = ifelse(str_detect(Replicate.Name, "atrix"), "Matrix", "Water"), sep = "_")) %>%
  select(Type,  MassFeature_Column, Area) %>%
  mutate(Area = as.numeric(as.character(Area))) %>%
  mutate(Area = ifelse(is.na(Area), 0, Area))%>%
  group_by(Type, MassFeature_Column) %>%
  summarise(Area = mean(Area, na.rm = TRUE)) %>% ungroup() %>%
  spread(., Type, Area) %>%
  as.data.frame() %>%
  mutate(RFratio = (Standards_Matrix-Water_Matrix)/Standards_Water) %>%
  select(MassFeature_Column, RFratio)

RFratio.RP <- dat %>%
  filter(Column == "RP",
         str_detect(Replicate.Name, "InH2O"))%>%
  mutate(RFratio = QE.RF.ratio) %>%
  select(MassFeature_Column, RFratio) %>%
  unique()
  
RFratio.TMAP <- dat %>%
  filter(MassFeature_Column == "I132.102R9.3_HILIC_X_HILIC") %>%
  mutate(RFratio = 1) %>%
  select(MassFeature_Column, RFratio) %>% unique()

RFratio.all <- rbind(RFratio.HILIC, RFratio.RP) %>% rbind(RFratio.TMAP)
write_csv(RFratio.all, "Intermediates/RFratio_all.csv")

rm(RFratio.HILIC, RFratio.RP, RFratio.TMAP)

#Make new RFs for compounds we didn't have standards for earlier dates using relative response factors----
#Grab all the uncorrected RF
dat.RF.bydate.uncorrected <- dat.RF.bydate

rf.ratio.matcher <- read_csv(rf.ratio.matcher.file)
dat.RF.bydate.matchers <- dat.RF.bydate.uncorrected %>% ungroup() %>%
  filter(Identification %in% rf.ratio.matcher$Identification_Matched) %>%
  rename(Identification_Matched = Identification) %>%
  rename(RF_matcher = RF) %>%
  select(Identification_Matched, date, RF_matcher) 

dat.RF.bydate.tofix.goodDate <- dat.RF.bydate.uncorrected %>%
  filter(Identification %in% rf.ratio.matcher$Identification) %>%
  filter(date %in% c("180702", "180716")) %>%
  left_join(rf.ratio.matcher %>% select(Column:Identification, Identification_Matched)) %>%
  left_join(dat.RF.bydate.matchers) %>%
  mutate(relativeRF = RF/RF_matcher)

relRFs <- dat.RF.bydate.tofix.goodDate %>% ungroup() %>%
  select(Identification, Identification_Matched, relativeRF)

dat.RF.bydate.tofix.badDate <- dat.RF.bydate.uncorrected %>% ungroup() %>%
  filter(Identification %in% rf.ratio.matcher$Identification) %>%
  filter(!date %in% c("180702", "180716")) %>%
  left_join(rf.ratio.matcher %>% select(Column:Identification, Identification_Matched)) %>%
  left_join(dat.RF.bydate.matchers) %>%
  left_join(relRFs) %>%
  mutate(Adjusted_RF = relativeRF*RF_matcher) %>%
  select(Identification, date, Adjusted_RF)
  
#Make a corrected dat.RF.bydate.uncorrected
dat.RF.bydate.corrected <- dat.RF.bydate.uncorrected %>%
  left_join(dat.RF.bydate.tofix.badDate) %>%
  mutate(RFFlag = ifelse(is.na(Adjusted_RF), NA, "used relative RF")) %>%
  mutate(RF = ifelse(is.na(Adjusted_RF), RF, Adjusted_RF)) %>%
  select(-Adjusted_RF)

write_csv(dat.RF.bydate.corrected, "Intermediates/RF_byDate_corrected.csv")

#Add it back to alldat-----
sampKey <- read_csv(samp.key.file) %>%
  select(RP_InjectionCorrection, Sample.Name)

#MF_info <- read.csv(id.file) %>%
#  filter(Confidence == 1) %>%
#  gather(sample, area, KM1513.125m_A:U9_B)
dat.RF.bydate.corrected.short <-  dat.RF.bydate.corrected %>% ungroup() %>%
  mutate(Date = as.character(date)) %>%
  select(MassFeature_Column, Date, RF, RFFlag)

MF.info.short <- MF_info %>% select(MassFeature_Column, Identification) %>% unique() %>%
  mutate(MassFeature_Column = MassFeature_Column %>% 
           str_replace_all("HILICNeg", "HILIC") %>%
           str_replace_all("HILICPos", "HILIC"))

datsmp.quan <- read_csv(dat.file) %>%
  separate(Run.Cmpd, sep = " ", into = c("Sample.Name")) %>%
  mutate(MF_Frac = MassFeature)%>%
  mutate(Date = as.character(Date)) %>%
  mutate(MassFeature_Column = paste(MF_Frac, Column, sep = "_X_")) %>%
  mutate(MassFeature_Column = MassFeature_Column %>% 
           str_replace_all("HILICNeg", "HILIC") %>%
           str_replace_all("HILICPos", "HILIC")) %>%
  filter(MassFeature_Column %in% RFratio.all$MassFeature_Column) %>%
  left_join(sampKey)%>%
  left_join(RFratio.all) %>%
  left_join(dat.RF.bydate.corrected.short) %>%
  left_join(MF.info.short) %>%
  mutate(umolinvial = Adjusted_Area/RF/RFratio*RP_InjectionCorrection)
  
#Get better data for compounds with matched internal standards----
IS_key <- read_csv(is.names.file) %>%
  rename(Identification = Matched_Compounds,
         IS = Internal_Standards)

dat1 <- read.csv(dat.hilicpos.file.nostds) %>%
  mutate(Column = "HILICPos") %>% select(-Protein)
dat2 <- read.csv(dat.hilicneg.file.nostds)  %>%
  mutate(Column = "HILICNeg")
dat3 <- read.csv(dat.cyano.file.nostds)  %>%
  mutate(Column = "RP")

IS_dat <- rbind(dat1, dat2) %>%
  rbind(dat3) %>%
  filter(as.character(`Precursor.Ion.Name`) %in% IS_key$IS) %>%
  mutate(IS_Area = Area,
         IS = `Precursor.Ion.Name`) %>%
  select(IS_Area, IS,  Replicate.Name) %>%
  left_join(IS_key) %>%
  rename(Sample.Name = Replicate.Name)

datsmp.quan.newIS <- datsmp.quan %>%
  filter(Identification %in% IS_key$Identification) %>%
  left_join(IS_dat) %>%
  mutate(umolinvial_IS = as.numeric(as.character(Area))
         /as.numeric(as.character(IS_Area))*Concentration_nM/1000) 

datsmp.quan.IScorrected <- datsmp.quan %>%
  filter(!Identification %in% IS_key$Identification) %>%
  bind_rows(datsmp.quan.newIS) %>%
  mutate(umolinvial = ifelse(is.na(umolinvial_IS), umolinvial, umolinvial_IS))%>%
  mutate(QuanFlag = ifelse(is.na(umolinvial_IS), NA, "Matched IS")) %>%
  select(MassFeature:umolinvial, QuanFlag)


#Calculate the concentration in the environment, not just in the vial -----
MetaDat <- read_csv(meta.data.file)

#Add in dilution factor from standard.info
quanDat2 <- datsmp.quan.IScorrected %>%
  left_join(MetaDat %>% select(SampID, Dilution_Factor, Volume, PC_ave, PN_ave)) %>%
  mutate(nmolinEnviroave = umolinvial*10^-6*400/Volume*1000*Dilution_Factor) %>%
  left_join(MF_info2 %>% select(Identification, Emperical.Formula)) %>%
  unique()

#Okay dokay, go get how many Carbons and Nitrogens there are here.------
quanDat3 <- quanDat2  %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"), 
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1, 
                    str_extract(Emperical.Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(nmolCave = nmolinEnviroave*C,
         nmolNave = nmolinEnviroave*N ) %>%
  mutate(percentCave = nmolCave/(PC_ave*1000)*100,
         percentNave = nmolNave/(PN_ave*1000)*100) %>%
  select(Identification, MassFeature_Column, Rank, RankPercent, SampID, nmolinEnviroave:nmolinEnviroave, nmolCave:percentNave)


#Cacluate mole fractions of each compound -------
TotalMoles <- quanDat3  %>%
  select(SampID, nmolCave, nmolNave, Identification) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM = sum(as.numeric(nmolCave), na.rm = TRUE),
            totalNmeasured_nM = sum(as.numeric(nmolNave), na.rm = TRUE))

quanDat4 <- quanDat3 %>%
  left_join(TotalMoles) %>%
  mutate(molFractionC = nmolCave/totalCmeasured_nM, 
         molFractionN = nmolNave/totalNmeasured_nM)

quanDatSum <- quanDat4 %>%
  group_by(MassFeature_Column, Identification) %>%
  summarise(nmolEnviromed = median(nmolinEnviroave, na.rm  = T),
            nmolEnviromin = min(nmolinEnviroave, na.rm  = T),
            nmolEnviromax = max(nmolinEnviroave, na.rm  = T),
            nmolCmed = median(nmolCave, na.rm  = T),
            nmolCmin = min(nmolCave, na.rm  = T),
            nmolCmax = max(nmolCave, na.rm  = T), 
            percentCmed = median(percentCave, na.rm = T), 
            percentCmin = min(percentCave, na.rm = T),  
            percentCmax =  max(percentCave, na.rm = T), 
            molFractionmed = median(molFractionC, na.rm = T),
            molFractionmin = min(molFractionC, na.rm = T),
            molFractionmax = max(molFractionC, na.rm = T)) %>%
  arrange(desc(molFractionmed))

quanDatWide <- quanDat4 %>%
  select(Identification, SampID, molFractionC) %>%
  spread(data = ., value = molFractionC, key = SampID)



write.csv(quanDatSum, "Intermediates/Quantified_MFSummary.csv")
write.csv(quanDat4, "Intermediates/Quantified_LongDat.csv")


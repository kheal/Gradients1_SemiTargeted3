library(tidyverse)
library(RCurl)

#Set your datafiles----
stds.file1 <- "RawOutput/Stds/HILIC_Pos_Stds_Normal_Skyline_withAOA.csv"
stds.file3 <- "RawOutput/Stds/HILIC_NEG_Stds_with_AOA.csv"
stds.file4 <- "RawOutput/Stds/RP_Stds_with_AOA.csv"
mf.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
stds.info <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"
RF.matcher.file <- "MetaData/RFs_relativeRFmatcher.csv"

#Read in and tidy up your stds data-----
#Connect the standards runs
stds.dat.HILICPos <- read_csv(stds.file1) %>%
  mutate(Column = "HILIC",
         z = 1) %>%
  mutate(`Precursor Ion Name` = `Precursor Ion Name` %>% str_replace(., "Beta_GlutamicAcid", "beta-Glutamic acid"))
  
stds.dat.HILICNeg <- read_csv(stds.file3) %>%
  mutate(Column = "HILIC",
         z = -1) 

stds.dat.RP <- read_csv(stds.file4) %>%
  mutate(Column = "RP",
         z = 1) %>%
  mutate(`Replicate Name` = str_replace(`Replicate Name`, "1800302", "180302"))
  
stds.dat.combo <- rbind(stds.dat.HILICPos, stds.dat.HILICNeg)%>% 
  rbind(stds.dat.RP) %>%
  mutate(Area = as.numeric(Area))

#Get RFs of all compounds----------
mf.dat <- read_csv(mf.file) %>% select(MassFeature_Column, Identification, Column, z, mz)

stds.dat.combo2 <- stds.dat.combo %>%
  rename(Identification = `Precursor Ion Name`,
         RunID = `Replicate Name`) %>%
  left_join(mf.dat, by = c("Identification", "Column", "z")) %>%
  select(MassFeature_Column, Identification, RunID, Area, Column, z) 

#Get dates, dump matrix samples, get mixes if applicable
stds.dat.combo3 <- stds.dat.combo2 %>%
  mutate(Date = str_extract(RunID, "^\\d{6}"))%>%
  filter(!str_detect(RunID, "Matrix")) %>%
  mutate(Mix = str_extract(RunID, "Mix\\d"))

#Get concentrations in order to get RFs
stds.basic.info <- read.csv(text = getURL(stds.info), header = T) %>%
  rename(Identification = Compound.Name,
         Concentration = `Conc..uM`) %>%
  select(Identification, Concentration, Column, z, Date.added, HILICMix) 

stds.dat.combo4 <- left_join(stds.dat.combo3, stds.basic.info, by = c("Identification", "Column", "z")) %>%
  filter(HILICMix == Mix | is.na(Mix)) %>% select(-Mix, -HILICMix)

#Get RFs-----
stds.dat.combo5 <- stds.dat.combo4 %>%
  mutate(RF = ifelse(Date.added < Date, as.numeric(Area)/Concentration, NA))

RFs.combo <- stds.dat.combo5 %>%
  group_by(Identification, MassFeature_Column, Date) %>%
  summarise(RFmax = max(RF),
            RFmin = min(RF), 
            RF = mean(RF, na.rm = TRUE))

#Get relative RFs of HILICPos compounds without injected standards----
#First calculate relative RFs for the compounds we need to for each date
rf.ratio.matcher <- read_csv(RF.matcher.file)

dat.RF.bydate.matchers <- RFs.combo %>% ungroup() %>%
  filter(Identification %in% rf.ratio.matcher$Identification_Matched) %>%
  rename(Identification_Matched = Identification) %>%
  rename(RF_matcher = RF) %>%
  select(Identification_Matched, Date, RF_matcher) 

dat.RF.goodDate.RFrelcalc <- RFs.combo %>%
  filter(Identification %in% rf.ratio.matcher$Identification) %>%
  filter(!is.na(RFmax))%>%
  left_join(rf.ratio.matcher %>% select(Identification, Identification_Matched), by = "Identification") %>%
  left_join(dat.RF.bydate.matchers, by = c("Date", "Identification_Matched")) %>%
  mutate(relativeRF = RF/RF_matcher)

relRFs <- dat.RF.goodDate.RFrelcalc %>% ungroup() %>%
  group_by(Identification, Identification_Matched) %>%
  summarise(relativeRFmax = max(relativeRF, na.rm = TRUE),
            relativeRFmin = min(relativeRF, na.rm = TRUE),
            relativeRF = mean(relativeRF, na.rm = TRUE))

dat.RF.bydate.tofix.badDate <- RFs.combo %>% ungroup() %>%
  filter(Identification %in% rf.ratio.matcher$Identification) %>%
  filter(is.na(RFmax))%>%
  left_join(rf.ratio.matcher %>% select(Identification, Identification_Matched), by = "Identification") %>%
  left_join(dat.RF.bydate.matchers, by = c("Date", "Identification_Matched"))%>%
  left_join(relRFs, by = c("Identification", "Identification_Matched")) %>%
  mutate(Adjusted_RF = relativeRF*RF_matcher,
         Adjusted_RFmax = relativeRFmax*RF_matcher, 
         Adjusted_RFmin = relativeRFmin*RF_matcher) %>%
  select(Identification, Date, Adjusted_RF, Adjusted_RFmax, Adjusted_RFmin)

RFs.combo.withRels <- RFs.combo %>%
  left_join(dat.RF.bydate.tofix.badDate, by = c("Identification", "Date"))%>%
  mutate(RFFlag = ifelse(is.na(Adjusted_RF), NA, "used relative RF")) %>%
  mutate(RF = ifelse(is.na(Adjusted_RF), RF, Adjusted_RF),
         RFmax = ifelse(is.na(Adjusted_RF), RFmax, Adjusted_RFmax),
         RFmin = ifelse(is.na(Adjusted_RF), RFmin, Adjusted_RFmin)) %>%
  select(-Adjusted_RF, -Adjusted_RFmax, -Adjusted_RFmin)

#Calculate the RF ratios when possible
stds.dat.rfratio <- stds.dat.combo %>%
  rename(Identification = `Precursor Ion Name`,
         RunID = `Replicate Name`) %>%
  left_join(mf.dat, by = c("Identification", "Column", "z")) %>%
  filter(str_detect(RunID, "^19")) %>%
  mutate(Mix = str_extract(RunID, "Mix\\d")) %>% 
  left_join(., stds.basic.info, by = c("Identification", "Column", "z")) %>%
  filter(HILICMix == Mix | is.na(Mix)) %>% select(-Mix, -HILICMix) 

stds.dat.rfratio2 <- stds.dat.rfratio %>%
  select(Identification, Area, RunID) %>%
  mutate(RunNumber = str_extract(RunID, "_\\d$")) %>%
  mutate(Date = str_extract(RunID, "^\\d{6}"))%>%
  mutate(RunType = ifelse(str_detect(RunID, "StdsMix\\dInH2O")|
                            str_detect(RunID, "StdsInH2O"), "Std_in_h2O", 
                          ifelse(str_detect(RunID, "StdsMix\\dInMatrix") |
                                   str_detect(RunID, "StdsInMatrix"), "Std_in_matrix",
                          "Matrix_in_h2O"))) %>%
  select(-RunID)%>%
  spread(key = RunType, value = Area )
  
stds.dat.rfratio3 <- stds.dat.rfratio2 %>%
  mutate(RFratio = (Std_in_matrix - Matrix_in_h2O)/ Std_in_h2O) %>%
  group_by(Identification) %>%
  summarise(RFratio = mean(RFratio))

#Put them together!
RFs.combo.withRels.withRFratio <- RFs.combo.withRels %>%
  left_join(stds.dat.rfratio3, by = "Identification") 


#Write out your results-----
write_csv(RFs.combo.withRels.withRFratio, "Intermediates/RFsandRFratios.csv")

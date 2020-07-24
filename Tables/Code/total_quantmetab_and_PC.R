#TO DO: calculate errors for MGL and KM samples
library(tidyverse)
library(fuzzyjoin)

#Read in your dat files
dat.filename <- "Intermediates/Quantified_LongDat_Enviro.csv"
PC.dat.filename <- "MetaData/PCPN/KOK1606_PCPN_UW_Preliminary_OSU_KRH.csv"
meta.dat.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"
  
#Check out data
dat <- read_csv(dat.filename)
dat.pc <- read_csv(PC.dat.filename)
meta.dat <- read_csv(meta.dat.filename)%>% select(SampID, Station_1, Cruise, latitude, longitude, Depth)
meta.dat.wPC <- read_csv(meta.dat.filename)%>% 
  select(Station_1, Cruise, latitude, longitude, Depth, PC_ave, PC_Stdev, PN_ave, PN_Stdev) %>% unique()

#Clean it up
dat.clean <- dat %>%
  left_join(meta.dat, by = "SampID") %>%
  mutate(Station_1 = ifelse(str_detect(SampID, "KM1513"), "KM1513", Station_1)) %>%
  mutate(Station_1 = ifelse(str_detect(SampID, "S7E41"), "KOK1606", Station_1)) %>%
  group_by(Station_1, Depth, Cruise) %>%
  summarise(nmolave.metabC = mean(totalCmeasured_nM),
            nmolstd.metabC = sd(totalCmeasured_nM), 
            nmolave.metabN = mean(totalNmeasured_nM),
            nmolstd.metabN = sd(totalNmeasured_nM), 
            latitude = mean(latitude))%>%
  mutate(nmolstd.metabC = ifelse(nmolstd.metabC== 0, NA, nmolstd.metabC)) %>%
  mutate(nmolstd.metabN = ifelse(nmolstd.metabN==0, NA, nmolstd.metabN))


#
PC.dat.summary <- dat.pc %>%
  mutate(latitude.round = round(Latitude, digits = 0)) %>%
  group_by(latitude.round) %>%
  summarise(nmolave.PC = mean(PC)*1000,
            nmolstd.PC = sd(PC)*1000, 
            nmolave.PN = mean(PN)*1000,
            nmolstd.PN = sd(PN)*1000, 
            latitude = mean(Latitude))

mesh.KOK.PC <- dat.clean %>%
  filter(Cruise == "KOK1606") %>%
  difference_left_join(PC.dat.summary, by = "latitude", max_dist = 0.5,
                       distance_col = "lat.difference") %>%
  filter(Station_1 != 7 | lat.difference < 0.47) %>%
  filter(Station_1 != 8 | lat.difference < 0.46) 

mesh.nonKOK.PC <- dat.clean %>%
  filter(Cruise != "KOK1606") %>%
  left_join(meta.dat.wPC, by = c("Depth", "Cruise"))%>%
  mutate(nmolave.PC = mean(PC_ave)*1000,
            nmolstd.PC = sd(PC_Stdev)*1000, 
            nmolave.PN = mean(PN_ave)*1000,
            nmolstd.PN = sd(PN_Stdev)*1000)

mesh.KOK.and.nonKOK <- mesh.KOK.PC %>%
  bind_rows(mesh.nonKOK.PC) %>%
  mutate(fraction.PC = nmolave.metabC/nmolave.PC) %>%
  mutate(fraction.PC.sd = ifelse(is.na(nmolstd.metabC), 
                                 fraction.PC*(nmolstd.PC/nmolave.PC),
                                 ifelse(is.na(nmolstd.PC),
                                        fraction.PC*(nmolstd.metabC/nmolave.PC),
                                 fraction.PC*(nmolstd.metabC/nmolave.PC+nmolstd.PC/nmolave.PC)))) %>%
  mutate(fraction.PN = nmolave.metabN/nmolave.PN) %>%
  mutate(fraction.PN.sd = ifelse(is.na(nmolstd.metabN), 
                                 fraction.PN*(nmolstd.PN/nmolave.PN),
                                 ifelse(is.na(nmolstd.PN),
                                        fraction.PC*(nmolstd.metabN/nmolave.PN),
                                 fraction.PN*(nmolstd.metabN/nmolave.PN+nmolstd.PN/nmolave.PN)))) 


mesh.KOK.PC.new.names <- mesh.KOK.and.nonKOK %>%
  ungroup() %>%
  mutate(nmolave.metabC.withstd = ifelse(is.na(nmolstd.metabC), as.character(round(nmolave.metabC, digits = 1)) , 
                                         paste0(round(nmolave.metabC, digits = 1), " (", round(nmolstd.metabC, digits = 1), ")"))) %>%
  mutate(nmolave.PC.withstd = ifelse(is.na(nmolave.PC), NA,
                                     ifelse(is.na(nmolstd.PC),  paste0(round(nmolave.PC, digits = -1)),
                                     paste0(round(nmolave.PC, digits = -1), " (", round(nmolstd.PC, digits = -1), ")")))) %>%
  mutate(percent.PC.withstd = ifelse(is.na(nmolave.PC), NA,
                                     ifelse(is.na(nmolstd.PC), paste0(round(fraction.PC*100, digits = 1)),
                                     paste0(round(fraction.PC*100, digits = 1), " (", round(fraction.PC.sd*100, digits = 1), ")")))) %>%
  mutate(nmolave.metabN.withstd = ifelse(is.na(nmolstd.metabN), as.character(round(nmolave.metabN, digits = 1)) , 
                                         paste0(round(nmolave.metabN, digits = 1), " (", round(nmolstd.metabN, digits = 2), ")"))) %>%
  mutate(nmolave.PN.withstd = ifelse(is.na(nmolave.PN), NA,
                                     ifelse(is.na(nmolstd.PN),paste0(round(nmolave.PN, digits = -1)),
                                     paste0(round(nmolave.PN, digits = -1), " (", round(nmolstd.PN, digits = -1), ")")) ))%>%
  mutate(percent.PN.withstd = ifelse(is.na(nmolave.PN), NA,
                                     ifelse(is.na(nmolstd.PN),paste0(round(fraction.PN*100, digits = 1)),
                                     paste0(round(fraction.PN*100, digits = 1), " (", round(fraction.PN.sd*100, digits = 2), ")")))) %>%
  select(Cruise, latitude.x, Depth,  
         nmolave.metabC.withstd, nmolave.PC.withstd, percent.PC.withstd, nmolave.metabN.withstd, nmolave.PN.withstd, percent.PN.withstd) %>%
  rename(`Total measured meatbolites (nmol C L-1)` = nmolave.metabC.withstd, 
         `Total measured meatbolites (nmol N L-1)` = nmolave.metabN.withstd,
         `Total measured particulates (nmol C L-1)` = nmolave.PC.withstd, 
         `Total measured particulates (nmol N L-1)` = nmolave.PN.withstd,
         `Metabolites as a portion of particulate matter (C, %)` = percent.PC.withstd,
         `Metabolites as a portion of particulate matter (N, %)` = percent.PN.withstd,
         Latitude = latitude.x, 
         `Depth (m)` = Depth)
  
write_csv( mesh.KOK.PC.new.names, "Tables/Manuscript_tables/SuppTables/Total_QuanMetabandPC.csv")

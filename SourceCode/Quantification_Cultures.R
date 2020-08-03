library(tidyverse)
library(RCurl)

#TO DO: maybe replace the ones with have with internal standards, would take a lot more work
#TO DO: make summaries like the old one so we can put it on the summary box plot
#TO DO: add in carbon estimates per cell so we can get %
#   mutate(nmolmetab_perC = intracell_conc_umolCL*BioVol_perFilter_uL/nmolC_filtered_final, by = "CultureID") 


#Set your datafiles
dat.file1 <- "Intermediates/Culture_Intermediates/combined_long_cultures.csv" 
mf.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
RF.file <- "Intermediates/RFsandRFratios.csv"
meta.file <- "MetaData/CultureMetaData.csv"
std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"

#Get good MFs
mfs.good <- read_csv(mf.file) %>% select(MassFeature_Column, Identification)

#Get RFs
RFs <- read_csv(RF.file)

#Read in your culture data, cull it so its just the compounds we can quantify-----------
dat <- read_csv(dat.file1) %>%
  left_join(read_csv(mf.file) %>% select(MassFeature_Column, Column, z), by = "MassFeature_Column") %>%
  mutate(Vol_recon_uL = 400) %>%
  filter(MassFeature_Column %in% RFs$MassFeature_Column) 

#Attach information about the RP volume if needed
meta.dat <- read_csv(meta.file) %>% 
  rename(ID_rep = CultureID) %>% 
  select(ID_rep, RP_stdVolpersmpVol, `Cells filtered`)
meta.dat2 <- read_csv(meta.file) %>% 
  rename(ID_rep = CultureID) %>%  
  select(ID_rep, nmolC_filtered_final) 

dat <- dat %>%
  left_join(meta.dat, by = "ID_rep")

#Attach RF and RFratio information
dat.withRFs <- dat %>%
  left_join(RFs, by = c("Date", "MassFeature_Column"))

#Calculate intracellular concentration----------
dat.with.Quan <- dat.withRFs %>%
  mutate(intracell_conc_umolL = QC_area/RF*Vol_recon_uL/BioVol_perFilter_uL*1/RFratio) %>%
  mutate(intracell_conc_umolL = ifelse(Column == "RP", intracell_conc_umolL*RP_stdVolpersmpVol, intracell_conc_umolL))

#Get C and N for each compound ---------
IngallsStandards <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name) %>%
  select(Identification, Emperical.Formula) %>% unique()

dat.with.Quan2 <- dat.with.Quan %>%
  left_join(IngallsStandards, by = "Identification") %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"), 
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1, 
                    str_extract(Emperical.Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(intracell_conc_umolCL = intracell_conc_umolL*C,
         intracell_conc_umolNL = intracell_conc_umolL*N ) 
 # mutate(percentCave = nmolCave/(PC_ave*1000)*100, #Add something like this once we've got estimates of C
 #        percentNave = nmolNave/(PN_ave*1000)*100) %>% #Add something like this once we've got estimates of C
 # select(Identification, MassFeature_Column, Rank, RankPercent, SampID, nmolinEnviroave:nmolinEnviroave, nmolCave:percentNave)


#Cacluate mole fractions of each compound -------
TotalMoles <- dat.with.Quan2  %>%
  select(ID_rep, intracell_conc_umolCL, intracell_conc_umolNL, Identification) %>%
  group_by(ID_rep) %>%
  summarise(totalCmeasured_uM = sum(as.numeric(intracell_conc_umolCL), na.rm = TRUE),
            totalNmeasured_uM = sum(as.numeric(intracell_conc_umolNL), na.rm = TRUE))

dat.with.Quan3 <- dat.with.Quan2 %>%
  left_join(TotalMoles, by = "ID_rep") %>%
  mutate(molFractionC = intracell_conc_umolCL/totalCmeasured_uM, 
         molFractionN = intracell_conc_umolNL/totalNmeasured_uM)

dat.with.Quan4 <- dat.with.Quan3 %>% 
  left_join(meta.dat2, by = "ID_rep") %>%
  mutate(nmolmetab_perC = intracell_conc_umolCL*BioVol_perFilter_uL/nmolC_filtered_final) 

quanDatSum <- dat.with.Quan4 %>%
  group_by(MassFeature_Column, Identification) %>%
  summarise(umol.intracell.med = median(intracell_conc_umolL, na.rm  = T),
            umol.intracell.min = min(intracell_conc_umolL, na.rm  = T),
            umol.intracell.max = max(intracell_conc_umolL, na.rm  = T),
            umol.intracellC.med = median(intracell_conc_umolCL, na.rm  = T),
            umol.intracellC.min = min(intracell_conc_umolCL, na.rm  = T),
            umol.intracellC.max = max(intracell_conc_umolCL, na.rm  = T), 
            #    percentCmed = median(percentCave, na.rm = T),   #Add this back in eventually
            #    percentCmin = min(percentCave, na.rm = T),   #Add this back in eventually
            #    percentCmax =  max(percentCave, na.rm = T),  #Add this back in eventually
            molFractionmed = median(molFractionC, na.rm = T),
            molFractionmin = min(molFractionC, na.rm = T),
            molFractionmax = max(molFractionC, na.rm = T)) %>%
  arrange(desc(molFractionmed))

quanDatWide <- dat.with.Quan4 %>%
  select(Identification, ID_rep, molFractionC) %>%
  spread(data = ., value = molFractionC, key = ID_rep)

#Exploring for specific numbers
explore <- dat.with.Quan3 %>%
  filter(str_detect(Org_Type, "aptophyte")) %>%
  filter(str_detect(Identification, "Homarine"))

#Write it out :)------
write_csv(dat.with.Quan4, "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv")
write_csv(quanDatSum, "Intermediates/Culture_Intermediates/Quantified_MFSummary_cultures.csv")


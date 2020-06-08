# id.file <- "RawOutput//NoIsos_allFractions_VolNormed_BMISd_CVFiltered_xsetWIDE_wMS2_wIDs.csv"
# id.manual.file <- "RawOutput/MFs_match_wManual.csv"
# std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"
# dat.hilicpos.file <- "RawOutput/HILICPos_IntegrationsBigPeaksWTargeted2.csv"
# dat.hilicneg.file <- "RawOutput/HILICNeg_IntegrationsBigPeaksWTargeted.csv"
# dat.cyano.file <- "RawOutput/CyanoAq_IntegrationsBigPeaksWTargeted.csv"
# samp.key.file <- "MetaData/Sample_key.csv"
# dat.file <- "Intermediates/Longdata.csv"
# is.names.file <- "MetaData/InternalStandardNames.csv"
# meta.data.file <- "MetaData/SampInfo_wMetaData.csv"
# 


###FUNCTION CALL GOES HERE
quantify.metabs <- function(id.file, id.manual.file, std.url, dat.hilicpos.file, 
                            dat.hilicneg.file, dat.cyano.file, samp.key.file,
                            dat.file, is.names.file, meta.data.file){


library(RCurl)

#Get info on the standards that we can quantify -----
IDsetc <- read.csv(id.file) %>%
  select(MF_Frac, Confidence, BestMatch)
MFs_Fix <- read_csv(id.manual.file) %>%
  select(MF_Frac, Column, Manual_Annote, Manual_Confidence) %>%
  left_join(IDsetc) %>%
  mutate(BestMatch = ifelse(str_detect(MF_Frac, "^I\\d"), as.character(BestMatch), as.character(MF_Frac))) %>%
  mutate(Confidence = ifelse(str_detect(MF_Frac, "^I\\d"), Confidence, 1)) %>%
  mutate(BestMatch = ifelse(is.na(Manual_Annote), as.character(BestMatch), as.character(Manual_Annote))) %>%
  mutate(Confidence = ifelse(is.na(Manual_Annote), Confidence, Manual_Confidence)) 

MF_info <- MFs_Fix %>%
  filter(Confidence == 1) %>%
  mutate(MassFeature_Column = paste(MF_Frac, Column, sep = "_X_")) %>%
  select(MassFeature_Column, MF_Frac, Column, BestMatch)

#Find which mix its in, what concentration its in that mix
IngallsStandards <- read.csv(text = getURL(std.url), header = T) 

# Use 1 = RF ratio for MTA since it's so late in the chromatogram----
IngallsStandards <- IngallsStandards %>%
  mutate(QE.RF.ratio = ifelse(Compound.Name == "Methylthioadenosine" | Compound.Name == "Trimethylammonium Propionate (TMAP)"
                              , 1, as.numeric(as.character(QE.RF.ratio))))

IngallsStandards2 <- IngallsStandards %>%
  mutate(Column = ifelse(Column == "RP", "CyanoAq", ifelse(z == -1, "HILICNeg", "HILICPos")),
         BestMatch = Compound.Name) %>%
  select(BestMatch, Column, QE.RF.ratio, HILICMix, Conc..uM, Emperical.Formula)

MF_info2 <- MF_info %>%
  left_join(IngallsStandards2) %>%
  mutate(Conc..uM = ifelse(MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos", 4, Conc..uM))%>%
  filter(!is.na(Conc..uM))

dat1 <- read.csv(dat.hilicpos.file) %>%
  mutate(Column = "HILICPos") %>% select(-Protein)
dat2 <- read.csv(dat.hilicneg.file)  %>%
  mutate(Column = "HILICNeg")
dat3 <- read.csv(dat.cyano.file)  %>%
  mutate(Column = "CyanoAq")

# Pull out just the standards
dat <- rbind(dat1, dat2) %>%
  rbind(dat3) %>%
  mutate(MF_Frac = Precursor.Ion.Name) %>%
  mutate(MassFeature_Column = paste(MF_Frac, Column, sep = "_X_")) %>%
  filter(MassFeature_Column %in% MF_info2$MassFeature_Column) %>%
  filter(str_detect(Replicate.Name, "Std")) %>%
  left_join(MF_info2)

#Get RF ratio for HILIC compounds - Use the data from 180702----
RFratio <- dat %>%
  filter(!Column == "CyanoAq") %>%
  filter(str_detect(Replicate.Name, "180702|180709")) %>%
  filter(MassFeature_Column != "I132.102R9.3_HILICPos_X_HILICPos") %>%
  mutate(Type = paste(Standards = ifelse(str_detect(Replicate.Name, as.character(HILICMix)), "Standards", "Water"),     
                      Matrix = ifelse(str_detect(Replicate.Name, "atrix"), "Matrix", "Water"), sep = "_")) %>%
  select(Type,  MassFeature_Column, Area) %>%
  mutate(Area = as.numeric(as.character(Area))) %>%
  mutate(Area = ifelse(is.na(Area), 0, Area))%>%
  spread(., Type, Area) %>%
  as.data.frame() %>%
  mutate(RFratio = (Standards_Matrix-Water_Matrix)/Standards_Water)

RPdatRFs <- dat %>%
  filter(Column == "CyanoAq" | MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos",
         str_detect(Replicate.Name, "InH2O"))%>%
  mutate(RFratio = QE.RF.ratio)

#Get RFs for all compounds - Use the data from all stds in water----
RFs <- dat %>%
  filter(!str_detect(Replicate.Name, "atrix")) %>%
  mutate(Area = ifelse(MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos", 0, as.numeric(as.character(Area)))) %>%
  filter(!is.na(Area))%>%
  filter(!(str_detect(Replicate.Name, "180702|180709") &  !str_detect(Replicate.Name, as.character(HILICMix)))) %>% 
  left_join(. ,dat %>% 
              select(BestMatch, MassFeature_Column) %>%
              unique()) %>%
  left_join(MF_info2 %>% select(BestMatch, MassFeature_Column, Conc..uM)) %>%
  mutate(RF = as.numeric(as.character(Area))/as.numeric(Conc..uM)) %>%
  group_by(MassFeature_Column, BestMatch) %>%
  summarise(RFmax = max(RF),
            RFmin = min(RF)) %>%
  mutate(RFmin = ifelse(MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos", 615089984, RFmin),
         RFmax = ifelse(MassFeature_Column == "I132.102R9.3_HILICPos_X_HILICPos", 615089984, RFmax))


AllRFratiosandRFs <- RFratio %>% 
  select(RFratio, MassFeature_Column) %>%
  rbind(RPdatRFs %>% select(RFratio,  MassFeature_Column)) %>%
  unique()%>%
  inner_join(RFs %>% select(MassFeature_Column, RFmax, RFmin, BestMatch ))

repeatsToDump <- c("I150.0413R8.14_HILICNeg_X_HILICNeg", "Adenosyl Homocysteine_X_CyanoAq",  
                   "Guanosine_X_HILICNeg", "Inosine_X_HILICNeg", "Vitamin B3_X_CyanoAq", "Agmatine_X_CyanoAq",
                   "Allopurinol_X_CyanoAq", "Glutathione Disulfide_X_CyanoAq")



#Add it back to alldat----
sampKey <- read_csv(samp.key.file) %>%
  select(RP_InjectionCorrection, Sample.Name)

datsmp <- read_csv(dat.file) %>%
  separate(Run.Cmpd, sep = " ", into = c("Sample.Name")) %>%
  mutate(MF_Frac = MassFeature)%>%
  mutate(MassFeature_Column = paste(MF_Frac, Column, sep = "_X_")) %>%
  filter(MassFeature_Column %in% AllRFratiosandRFs$MassFeature_Column) %>%
  filter(!MassFeature_Column %in% repeatsToDump ) %>%
  left_join(sampKey) %>%
  left_join(AllRFratiosandRFs) %>%
  mutate(RFratio = as.numeric(RFratio))

quanDat <- datsmp %>%
  mutate(RFave = as.numeric(rowMeans(datsmp[, c("RFmin", "RFmax")]))) %>%
  mutate(umolinvialave = Adjusted_Area/RFave/RFratio*RP_InjectionCorrection,
         umolinvialmax = Adjusted_Area/RFmin/RFratio*RP_InjectionCorrection,
         umolinvialmin = Adjusted_Area/RFmax/RFratio*RP_InjectionCorrection) 

#Get better data for compounds with matched internal standards----
IS_key <- read_csv(is.names.file) %>%
  rename(MF_Frac = Matched_Compounds,
         IS = Internal_Standards)

dat1 <- read.csv(dat.hilicpos.file) %>%
  mutate(Column = "HILICPos") %>% select(-Protein)
dat2 <- read.csv(dat.hilicneg.file)  %>%
  mutate(Column = "HILICNeg")
dat3 <- read.csv(dat.cyano.file)  %>%
  mutate(Column = "CyanoAq")

IS_dat <- rbind(dat1, dat2) %>%
  rbind(dat3) %>%
  filter(as.character(`Precursor.Ion.Name`) %in% IS_key$IS) %>%
  mutate(IS_Area = Area,
         IS = `Precursor.Ion.Name`) %>%
  select(IS_Area, IS,  Replicate.Name) %>%
  left_join(IS_key) %>%
  rename(BestMatch = MF_Frac)

smp_dat <- rbind(dat1, dat2) %>%
  rbind(dat3) %>%
  mutate(MF_Frac = Precursor.Ion.Name)%>%
  filter(str_detect(Replicate.Name, "mp_"))%>%
  mutate(MassFeature_Column = paste(MF_Frac, Column, sep = "_X_")) %>%
  left_join(quanDat %>% select( MassFeature_Column, BestMatch) %>% unique()) %>%
  filter(BestMatch %in% IS_key$MF_Frac,
         !is.na(BestMatch))%>%
  left_join(IS_dat)%>%
  mutate(umolinvial_IS = as.numeric(as.character(Area))
         /as.numeric(as.character(IS_Area))*Concentration_nM/1000) %>%
  select(Replicate.Name, MF_Frac, umolinvial_IS, BestMatch) %>% unique() 

smp_dat_tobind <- smp_dat %>%
  mutate(SampID = str_replace_all(Replicate.Name, "^\\d\\d\\d\\d\\d\\d_Smp_", "")) %>%
  select(BestMatch, SampID, umolinvial_IS)


#Put them together-----
quanDat1 <- quanDat %>%
  left_join(smp_dat_tobind) %>%
  mutate(umolinvialave = ifelse(is.na(umolinvial_IS), umolinvialave, umolinvial_IS), 
         umolinvialmax = ifelse(is.na(umolinvial_IS), umolinvialmax, NA),
         umolinvialmin = ifelse(is.na(umolinvial_IS), umolinvialmin, NA))

#Get MetaData about each of the samples----
MetaDat <- read_csv(meta.data.file)


#Add in dilution factor from standard.info
quanDat2 <- quanDat1 %>%
  left_join(MetaDat %>% select(SampID, Dilution_Factor, Volume, PC_ave, PN_ave)) %>%
  mutate(nmolinEnviroave = umolinvialave*10^-6*400/Volume*1000*Dilution_Factor) %>%
  left_join(MF_info2 %>% select(BestMatch, Emperical.Formula)) %>%
  unique()

#Okay dokay, go get how many Carbons and Nitrogens there are here.
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
  select(BestMatch, MassFeature_Column, Rank, RankPercent, SampID, nmolinEnviroave:nmolinEnviroave, nmolCave:percentNave)

quanDatSum <- quanDat3 %>%
  group_by(MassFeature_Column, BestMatch) %>%
  summarise(nmolEnviromed = median(nmolinEnviroave, na.rm  = T),
            nmolEnviromin = min(nmolinEnviroave, na.rm  = T),
            nmolEnviromax = max(nmolinEnviroave, na.rm  = T),
            nmolCmed = median(nmolCave, na.rm  = T),
            nmolCmin = min(nmolCave, na.rm  = T),
            nmolCmax = max(nmolCave, na.rm  = T), 
            percentCmed = median(percentCave, na.rm = T), 
            percentCmin = min(percentCave, na.rm = T),  
            percentCmax =  max(percentCave, na.rm = T)) %>%
  arrange(desc(nmolEnviromed))


#Cacluate mole fractions of each compound -------
TotalMoles <- quanDat3  %>%
  select(SampID, nmolCave, nmolNave, BestMatch) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM = sum(as.numeric(nmolCave), na.rm = TRUE),
            totalNmeasured_nM = sum(as.numeric(nmolCave), na.rm = TRUE))

quanDat4 <- quanDat3 %>%
  left_join(TotalMoles) %>%
  mutate(molFractionC = nmolCave/totalCmeasured_nM, 
         molFractionN = nmolCave/totalNmeasured_nM)

quanDatSum <- quanDat4 %>%
  group_by(MassFeature_Column, BestMatch) %>%
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
  select(BestMatch, SampID, molFractionC) %>%
  spread(data = ., value = molFractionC, key = SampID)

datToReturn <- list(quanDat4, RFratio, RFs )

return(datToReturn)

}


#write.csv(quanDatSum, "Quantified_MFSummary.csv")
#write.csv(quanDat4, "Quantified_LongDat.csv")
#write.csv(quanDatWide, "Quantified_WideDat.csv")


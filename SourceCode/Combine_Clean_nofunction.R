library(tidyverse)
library(here)
library(RCurl)

#List of compounds I've seen in multiple fractions (SAM, adenine etc) or that were only observed in a few environmental samples or that are repeats of eachother 
repeats.to.dump <- c("Adenosyl Homocysteine_X_RP", 
                     "I116.0706R7.18_HILICNeg_X_HILICNeg", 
                     "Deoxyadenosine_X_RP", 
                     "I252.1093R1.15_CyanoAq_X_RP",
                     "I146.045R12.11_HILICNeg_X_HILICNeg", 
                     "I130.05R12.12_HILICPos_X_HILICPos",
                     "I128.0343R12.12_HILICNeg_X_HILICNeg", 
                     "I152.0567R0.9_CyanoAq_X_RP", 
                     "I150.0413R8.14_HILICNeg_X_HILICNeg", 
                     "I284.0992R1.08_CyanoAq_X_RP",
                     "Inosine_X_HILICPos", 
                     "Kynurenine_X_HILICPos",
                     "I132.1021R0.87_CyanoAq_X_RP", 
                     "I130.0864R6.66_HILICNeg_X_HILICNeg",
                     "I136.0395R7.72_HILICNeg_X_HILICNeg", 
                     "Vitamin B3_X_RP", 
                     "Guanosine_X_HILICNeg",
                     "I126.0222R10.87_HILICPos_X_HILICPos", 
                     "I206.9967R1.94_HILICNeg_X_HILICNeg",
                     "I284.2193R6.29_CyanoAq_X_RP", 
                     "I340.2597R6.69_CyanoAq_X_RP", 
                     "I567.4317R6.29_CyanoAq_X_RP", 
                     "I679.5123R6.69_CyanoAq_X_RP",
                     "I155.0815R6.1_HILICPos_X_HILICPos",
                     "I137.0459R1.2_CyanoAq_X_RP", 
                     "I136.0619R0.89_CyanoAq_X_RP", 
                     "I268.1042R1.02_CyanoAq_X_RP")


# Define all your inputs 
cyano.dat.file <- "Intermediates/BMISReports/BMISd_Areas_long_Cyano.csv"
hilic.dat.file <- "Intermediates/BMISReports/BMISd_Areas_long_HILIC.csv"
sampkey <- "MetaData/Sample_key.csv"
RSD_cut.off <- 0.5 #This is max RSD allowed in multiple injections of the pooled sample
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
dat.file.wide <- "Intermediates/WideArea.csv" 
id.manual.file <- "RawOutput/MFs_match_wManual.csv"
std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"

##Import  longData, turn it into wide data, wide ranked data, wide ranked percentile.-----
dat.cyano <- read.csv(cyano.dat.file) %>% mutate(Column = "RP")
dat.hilic <- read.csv(hilic.dat.file) 
  
dat <- rbind(dat.cyano, dat.hilic) %>%
  mutate(SampID = paste(SampID, replicate, sep = "_")) %>%
  filter(type == "Smp") %>%
  mutate(Area = ifelse(is.na(Area), 0, Area)) %>%
  mutate(MassFeature = MassFeature %>% str_replace(., "Beta_GlutamicAcid", "beta-Glutamic acid"))

meta.dat <- read_csv(meta.dat.file) %>% select(SampID, Volume)
  
#Dump bad MFs and bad samples here, change the name for beta-Glutamic acid name (from Beta_GlutamicAcid) here----
dat.dumped <- dat %>%
  mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
  filter(!str_detect(SampID, "S2C1" )) %>%
  filter(!MassFeature_Column %in% repeats.to.dump) 

dat.dumped <- dat.dumped %>% left_join(meta.dat) %>%
    mutate(Adjusted_Area_VolNormed = Adjusted_Area/Volume)
  
#Calculate rank and rank percent----
Columns <- unique(dat.dumped$Column)
column_list <- list(Columns)
  for (j in (1:length(Columns))){
    datsub <- dat.dumped %>%
      filter(Column == Columns[j])
    Samps <- unique(datsub$SampID)
    dat_rank_list <- list()
    for (i in (1:length(Samps))){
      dat_rank_list[[i]] <- datsub %>%
        filter(SampID == Samps[i]) %>%
        mutate(Rank = dense_rank(desc(Area)),
               RankPercent = (1-Rank/max(Rank)))}
    column_list[[j]] <- do.call(rbind, dat_rank_list)
  }
dat_rank_all <- do.call(rbind, column_list) 
write_csv(as.data.frame(dat_rank_all), "Intermediates/Longdata.csv")

#Make a wide df of the area----
datWide_Area <- dat_rank_all  %>%
    select(SampID,  MassFeature_Column, Adjusted_Area_VolNormed) %>%
    spread(., SampID, Adjusted_Area_VolNormed) %>%
    as.data.frame()
write_csv(datWide_Area, "Intermediates/WideArea.csv")

#Make a wide df of the rank----
datWide_Rank  <- dat_rank_all  %>%
    mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
    select(SampID, MassFeature_Column, Rank) %>%
    spread(., SampID, Rank) %>%
    as.data.frame()
write_csv(datWide_Rank, "Intermediates/WideRank.csv")

#Make a wide df of the rank percent----
datWide_RankPercent  <- dat_rank_all  %>%
    mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
    select(SampID, MassFeature_Column, RankPercent) %>%
    spread(., SampID, RankPercent) %>%
    as.data.frame()
write_csv(datWide_RankPercent, "Intermediates/WideRankPercentile.csv")

#Clean up the IDs and the names on the intermediate files--------


# dat.file <-"Intermediates/WideArea_ModID.csv" 
# id.manual.file <- "RawOutput/MFs_match_wManual.csv"
# std.url = "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"


dat <- read_csv(dat.file.wide)
mfs.info <- read_csv(id.manual.file) %>%
  mutate(Column = Column %>% str_replace_all("CyanoAq", "RP")) %>%
  mutate(MassFeature_Column = paste(MF_Frac, Column, sep = "_X_"))

mfs.to.fix.2 <-  mfs.info %>%
   select(MF_Frac, Manual_Annote, Manual_Confidence)  %>%
  filter(!is.na(Manual_Annote))

mfs.info.2 <-  mfs.info %>%
  select(MF_Frac, Column, BestMatch, Confidence)  %>%
  left_join(mfs.to.fix.2) 

dat.2 <- dat %>%
  select(MassFeature_Column) %>%
  mutate(MassFeature_Column_copy = MassFeature_Column) %>%
  separate(MassFeature_Column_copy, into = c("MF_Frac", "Column"), sep = "_X_") %>%
  left_join(mfs.info.2) %>%
  mutate(BestMatch = ifelse(is.na(Manual_Annote), BestMatch, Manual_Annote)) %>%
  mutate(Confidence = ifelse(is.na(Manual_Annote), Confidence, Manual_Confidence)) %>%
  mutate(BestMatch = ifelse(BestMatch == "unknown", NA, BestMatch) ) %>%
  mutate(Identification = ifelse(BestMatch == "Unknown", NA, BestMatch) ) 
  
dat.3 <- dat.2 %>%
  select(MassFeature_Column, Identification, Confidence) %>%
  left_join(mfs.info %>% select(MassFeature_Column, mz, rt, Column)) %>%
  mutate(z = ifelse(Column == "HILICNeg", -1, 1)) %>%
  mutate(Column = ifelse(Column == "RP", Column, "HILIC"))

dat.3.sub <- dat.3 %>%
  filter(is.na(mz))

Names <- dat.3.sub$Identification

std.dat.sub <- read.csv(text = getURL(std.url), header = T) %>%
  select(Compound.Name, m.z, RT..min., z, Column) %>%
  filter(Compound.Name %in% Names) %>%
  rename(Identification = Compound.Name)

dat.3.sub.2 <- dat.3.sub %>%
  left_join(std.dat.sub) %>%
  mutate(mz = m.z, rt = RT..min.*60) %>% select(-RT..min., -m.z)
  
dat.4 <-dat.3 %>%
  filter(!is.na(mz)) %>% bind_rows(dat.3.sub.2)

dat.final <- dat %>% left_join(dat.4) %>%
  select(colnames(dat.4), everything())

write_csv(dat.final, "Intermediates/WideArea_withIDinfo.csv")



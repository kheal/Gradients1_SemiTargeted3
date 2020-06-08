library(tidyverse)

#Set your datafiles

# Set your output files that have been QCd in Skyline
dat.file1 <- "Intermediates/Culture_Intermediates/QCd_HILICNeg.csv"
dat.file2 <- "Intermediates/Culture_Intermediates/QCd_HILICPos.csv"
dat.file3 <- "Intermediates/Culture_Intermediates/QCd_RP.csv"
meta.dat.file <- "MetaData/CultureMetaData.csv"

# Read in your QCd files, make appropriate mods
dat.HILICNeg <- read_csv(dat.file1, skip = 1) %>%
  mutate(MassFeature_Column = paste(Precursor.Ion.Name, "HILICNeg", sep = "_X_"))

dat.HILIPos <- read_csv(dat.file2, skip = 1) %>%
  mutate(MassFeature_Column = paste(Precursor.Ion.Name, "HILICPos", sep = "_X_"))

dat.RP <- read_csv(dat.file3, skip = 1)  %>%
  mutate(MassFeature_Column = paste(Precursor.Ion.Name, "RP", sep = "_X_"))

# Combine.  
dat1 <- dat.HILICNeg %>%
  rbind(dat.HILIPos) %>%
  rbind(dat.RP) %>%
  select(Replicate.Name, Precursor.Ion.Name, MassFeature_Column, QC_area) %>%
  mutate(Replicate.Name = Replicate.Name %>% str_replace_all(., "P5-5", "P55"))

#Meta Dat
meta.dat <- read_csv(meta.dat.file) %>% rename(ID_rep = CultureID) %>%
  select(ID_rep:Species)

# Seperate your data, get rid of internal standards, add in biovolume normalization
dat2 <- dat1 %>%
  separate(Replicate.Name, into = c("Date", "type", "ID", "replicate"), sep = "_", extra = "drop") %>%
  mutate(ID_rep = paste(ID, replicate, sep = "_")) %>%
  filter(!str_detect(Precursor.Ion.Name,"13C"),
         !str_detect(Precursor.Ion.Name,"15N"),
         !str_detect(Precursor.Ion.Name,"eavy"),
         !str_detect(Precursor.Ion.Name,"D-"),
         !str_detect(Precursor.Ion.Name,"d-IAA"),
         !str_detect(Precursor.Ion.Name,"d4"),
         !str_detect(Precursor.Ion.Name,"D5"),
         !str_detect(Precursor.Ion.Name,"D6-"),
         !str_detect(Precursor.Ion.Name,"D8"),
         !str_detect(Precursor.Ion.Name,"D3"),
         !str_detect(Precursor.Ion.Name,"D7"),
         !str_detect(Precursor.Ion.Name,"D4"),
         !str_detect(Precursor.Ion.Name, "d3-")) %>%
  mutate(PresAbs = ifelse(QC_area < 1 |is.na(QC_area), 0, 1)) %>%
  mutate(LogArea = ifelse(QC_area < 1 |is.na(QC_area), 0, log10(QC_area))) %>%
  select(-type) %>% left_join(meta.dat, by = "ID_rep") %>%
  mutate(biovolArea = QC_area/BioVol_perFilter_uL) %>%
  mutate(LogBioArea = ifelse(biovolArea < 0 |is.na(biovolArea), 0, log10(biovolArea)))


# Make wide dataframes for rawArea, presabs, and LogArea
dat.wide.rawArea <- dat2 %>%
  select(ID_rep, MassFeature_Column, QC_area) %>%
  spread(., ID_rep, QC_area) 

dat.wide.presabs <- dat2 %>%
  select(ID_rep, MassFeature_Column, PresAbs) %>%
  spread(., ID_rep, PresAbs) 

dat.wide.LogArea <- dat2 %>%
  select(ID_rep, MassFeature_Column, LogArea) %>%
  spread(., ID_rep, LogArea) 

dat.wide.biovolArea <- dat2 %>%
  select(ID_rep, MassFeature_Column, biovolArea) %>%
  spread(., ID_rep, biovolArea) 

dat.wide.LogBioArea <- dat2 %>%
  select(ID_rep, MassFeature_Column, LogBioArea) %>%
  spread(., ID_rep, LogBioArea) 


write_csv(dat2, "Intermediates/Culture_Intermediates/combined_long_cultures.csv")
write_csv(dat.wide.LogBioArea, "Intermediates/Culture_Intermediates/combined_wide_LogBioArea.csv")

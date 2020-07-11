library(here)
library(tidyverse)
library(xtable)
library(lubridate)

samp.info <- read_csv("MetaData/CultureMetaData.csv")
samp.info.short <- samp.info %>%
  select(CultureID, Species, strain, Org_Type, CultureID_short, BioVol_perFilter_uL, volume_reconst, HILIC_injectVol, RP_injectVolume, nmolC_filtered_final, Method) %>%
  rename(`Sample ID` = CultureID,
         `Broad taxon` = Org_Type, Strain = strain, 
         `Short ID` = CultureID_short, 
         `Biovolume filtered` = BioVol_perFilter_uL, 
         `Reconstitution volume` = volume_reconst,
         `injection volume (HILIC)` = HILIC_injectVol, 
         `injection volume (RP)` = RP_injectVolume, 
         `nmol C filtered` = nmolC_filtered_final,
         `C estimation method` = Method) %>%
  arrange(Species) %>%
  arrange(`Broad taxon`) 

write_csv(samp.info.short, "Tables/Manuscript_tables/SuppTables/Culture_descriptions_full.csv")

library(tidyverse)

#Read in your dat files
dat.filename <- "Intermediates/WideArea_withIDinfo.csv"


#Check out data
dat <- read_csv(dat.filename)
dat.names <- read_csv(dat.namematcher.filename)

#Clean it up
dat.clean <- dat %>%
  select(SampID,  cruise, latitude, longitude, Depth, UTC_real, SampType, Volume) %>%
  rename(`Sample ID` = SampID,
         `Cruise ID` = cruise,
         Latitude = latitude, 
         Longitude =  longitude,
         `Depth (m)` = Depth,
         `Time (UTC)` = UTC_real,
         `Sampling Method` = SampType,
         `Volume (L)` = Volume, 
          `Dilution factor` = Dilution_Factor) %>%
  filter(`Sample ID` %in% colnames(dat.names))

write_csv( dat.clean, "Tables/Manuscript_tables/SuppTables/Full_Enviro_Samp_Descriptions.csv")

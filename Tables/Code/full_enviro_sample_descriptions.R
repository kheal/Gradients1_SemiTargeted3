library(tidyverse)

#Read in your dat files
dat.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"
dat.namematcher.filename <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"

#Check out data
dat <- read_csv(dat.filename)
dat.names <- read_csv(dat.namematcher.filename)

#Get binned latitude
dat.lat.bin <- dat %>%
  select(SampID, Station_1, latitude) %>%
  filter(!is.na(Station_1)) %>%
  group_by(Station_1) %>%
  summarise(binned_lat = mean(latitude))

#Clean it up
dat.clean <- dat %>%
  left_join(dat.lat.bin) %>%
  mutate(SampType2 = ifelse(SampType == "S", "niskin", "underway"))%>%
  select(SampID,  cruise, latitude, longitude, binned_lat, Depth, UTC_real, 
         SampType2, Volume) %>%
  rename(`Sample ID` = SampID,
         `Cruise ID` = cruise,
         Latitude = latitude, 
         Longitude =  longitude,
         `Binned latitude` = binned_lat,
         `Depth (m)` = Depth,
         `Time (UTC)` = UTC_real,
         `Sampling Method` = SampType2,
         `Volume (L)` = Volume) %>%
  filter(`Sample ID` %in% colnames(dat.names))

write_csv( dat.clean, "Tables/Manuscript_tables/SuppTables/Full_Enviro_Samp_Descriptions.csv")

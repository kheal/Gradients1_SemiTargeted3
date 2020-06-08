# %TO DO: make supplmental Table of full description of samples
# \begin{table}[ht]
# \centering
# \caption{\label{FullEnviroSampleDescriptions} Full sample descriptions for environmental samples} 
# \end{table}

library(tidyverse)

#Read in your dat files
dat.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"
dat.namematcher.filename <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"

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
         `Volume (L)` = Volume) %>%
  filter(`Sample ID` %in% colnames(dat.names))

#Write out appropriate comment
comment <- "Sample descriptions for individual environmental samples collected and anlayzed"

con <- file("Tables/Manuscript_tables/SuppTables/Full_Enviro_Samp_Descriptions.csv", open="wt")
writeLines(paste(comment), con)
write.csv( dat.clean, con)
close(con)

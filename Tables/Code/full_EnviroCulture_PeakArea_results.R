library(tidyverse)

#Read in your dat files
dat.filename.enviro <- "Intermediates/WideArea_withIDinfo.csv"
dat.filename.cultures <- "Intermediates/Culture_Intermediates/combined_long_cultures.csv"


#Check out data
dat.enviro <- read_csv(dat.filename.enviro)
dat.cultures <- read_csv(dat.filename.cultures) %>%
  select(MassFeature_Column, ID_rep, biovolArea) %>%
  pivot_wider(id_cols = MassFeature_Column, names_from = ID_rep, values_from = biovolArea)
dat.combo <- dat.enviro %>% left_join(dat.cultures, by = "MassFeature_Column")

write_csv( dat.combo, "Tables/Manuscript_tables/SuppTables/Full_EnviroCulture_PeakAreas.csv")

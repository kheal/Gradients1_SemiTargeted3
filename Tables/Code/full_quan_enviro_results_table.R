#Add better names

library(tidyverse)
library(RCurl)

#Read in your dat files
dat.filename <- "Intermediates/Quantified_LongDat_Enviro.csv"

#Get better names for standards
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         `Compound name in figures` = Compound.Name_figure,
         `Complete compound name` = Compound.Name) %>%
  select(`Compound name in figures`, `Complete compound name`, Identification) %>% unique()

#Check out data
dat <- read_csv(dat.filename)


#Make it wide
dat.wide <- dat %>%
  select(Identification, SampID, nmolinEnviroave) %>%
  spread(key = SampID, value = nmolinEnviroave) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  select(`Compound name in figures`, `Complete compound name`, everything())

write_csv(dat.wide, "Tables/Manuscript_tables/SuppTables/Full_Quan_Enviro_Results.csv")

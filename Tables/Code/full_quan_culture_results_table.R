#Add better names

library(tidyverse)
library(RCurl)

#Read in your dat files
dat.filename <- "Intermediates/Culture_Intermediates/combined_long_withquan.csv"

#Get better names for standards
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         `Compound name in figures` = Compound.Name_figure,
         `Complete compound name` = Compound.Name) 

#Check out data
dat <- read_csv(dat.filename)


#Make it wide
dat.wide <- dat %>%
  unique() %>%
  filter(!is.na(Identification))%>%
  select(Identification, ID_rep, intracell_conc_umolL) %>%
  spread(key = ID_rep, value = intracell_conc_umolL) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  select(`Compound name in figures`, `Complete compound name`, everything())

write_csv(dat.wide, "Tables/Manuscript_tables/SuppTables/Full_Quan_Culture_Results.csv")

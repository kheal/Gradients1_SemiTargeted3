# %TO DO: make supplemental table of quantification scheme for each compound 
# \begin{table}[ht]
# \centering
# \caption{\label{Quantified_methods} Summary of quantification method used for each quantifiable compound.} 
# \end{table}


library(tidyverse)

#Read in your dat files
dat.filename <- "Intermediates/Quantified_LongDat_Enviro.csv"
dat.filename.culture <- "Intermediates/Culture_Intermediates/combined_long_withquan.csv"
dat.matcher.filename <- "MetaData/Date_sampleset_key.csv"
dat.RF.matcher.filename <- "MetaData/RFs_relativeRFmatcher.csv"

#Check out data
dat <- read_csv(dat.filename)
dat.cul <- read_csv(dat.filename.culture)
date.matcher <- read_csv(dat.matcher.filename)
RF.matcher <- read_csv(dat.RF.matcher.filename) 

#Make Enviro data into good format
dat.enivro <- dat %>%
  select(Identification, RFFlag, Date) %>%
  left_join(date.matcher %>% select(Date, Sample_set), by = "Date") %>%
  select(-Date) %>% unique() %>%
  arrange(Identification) %>%
  mutate(Quantification_method = ifelse(RFFlag == "used relative RF", "relative RF and RF ratio", 
                                         ifelse(RFFlag == "Matched IS", "isotopologue", NA) )) %>%
  mutate(Quantification_method = ifelse(is.na(Quantification_method), "external RF and RF ratio",  Quantification_method)) %>%
  left_join(RF.matcher, by = "Identification") %>%
  mutate(`Proxy compound for relative RF` = ifelse(Quantification_method == "relative RF and RF ratio", 
                                                   Identification_Matched, NA)) %>%
  select(Identification, Quantification_method, `Proxy compound for relative RF`, Sample_set)

dat.enivro.clean <- dat.enivro %>%
  group_by(Identification, Quantification_method, `Proxy compound for relative RF`) %>%
  summarise(Sample_set_Enviro = paste(Sample_set, collapse=", ")) %>%
  filter(!is.na(Identification))

#Make Culture data into good format
dat.culture <- dat.cul %>%
  select(Identification, RFFlag, Date) %>%
  left_join(date.matcher %>% select(Date, Sample_set), by = "Date")%>%
  select(-Date) %>% unique() %>%
  arrange(Identification) %>%
  mutate(Quantification_method = ifelse(RFFlag == "used relative RF", "relative RF and RF ratio", 
                                        ifelse(RFFlag == "Matched IS", "isotopologue", NA) )) %>%
  mutate(Quantification_method = ifelse(is.na(Quantification_method), "external RF and RF ratio",  Quantification_method)) %>%
  left_join(RF.matcher, by = "Identification") %>%
  mutate(`Proxy compound for relative RF` = ifelse(Quantification_method == "relative RF and RF ratio", 
                                                   Identification_Matched, NA)) 

dat.culture.clean <- dat.culture %>%
  group_by(Identification, Quantification_method, `Proxy compound for relative RF`) %>%
  summarise(Sample_set_Cultures = paste(Sample_set, collapse=", ")) %>%
  filter(!is.na(Identification))


#Clean it up!
dat.clean <- dat.enivro.clean %>%
  full_join(dat.culture.clean,
            by = c("Identification", "Quantification_method", "Proxy compound for relative RF")) %>%
  rename(Compound = Identification,
         `Quantification method` = Quantification_method,
         `Environmental samples applied` = Sample_set_Enviro, 
         `Culture samples applied` = Sample_set_Cultures) %>%
  arrange(Compound)

#Write out appropriate comment
comment <- "Quantification method for each quantified metabolite in each sample set.  Proxy compound (when applicable) is the compound by which a relative response factor (RF) was calculated"

con <- file("Tables/Manuscript_tables/SuppTables/Full_Quan_Methods.csv", open="wt")
writeLines(paste(comment), con)
write.csv(dat.clean, con)
close(con)

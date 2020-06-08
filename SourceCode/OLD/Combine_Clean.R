# Define all your inputs 
# cyano.dat.file = "Intermediates/BMISd_Areas_long_Cyano.csv"
# hilic.dat.file = "Intermediates/BMISd_Areas_long_HILIC.csv"
# sampkey <- "MetaData/Sample_key.csv"
# RSD_cut.off <- 0.5 #This is max RSD allowed in multiple injections of the pooled sample
# meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"

#Things to return 
#Normed_Samp_datLong is a data frame with bionormed data of just the samples

Combine_and_Clean <- function(cyano.dat.file, hilic.dat.file, repeats.to.dump, meta.dat.file ){

  library(tidyverse)
  CombineReport <-list()
  ##Import  longData, turn it into wide data, wide ranked data, wide ranked percentile.
  dat1 <- read.csv(cyano.dat.file) %>% mutate(Column = "RP")
  dat2 <- read.csv(hilic.dat.file) 
  
  dat <- rbind(dat1, dat2) %>%
    mutate(SampID = paste(SampID, replicate, sep = "_")) %>%
    filter(type == "Smp") %>%
    mutate(Area = ifelse(is.na(Area), 0, Area))
  meta.dat <- read_csv(meta.dat.file) %>% select(SampID, Volume)
  
  dat <- dat %>% left_join(meta.dat) %>%
    mutate(Adjusted_Area_VolNormed = Adjusted_Area/Volume)
  
  Columns <- unique(dat$Column)
  column_list <- list(Columns)
  for (j in (1:length(Columns))){
    datsub <- dat %>%
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
  dat_rank_all <- do.call(rbind, column_list) %>%
    mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
    filter(!MassFeature_Column %in% repeats.to.dump) %>%
    filter(!str_detect(SampID, "S2C1" ))
  
  
######  write_csv(dat_rank_all, "Longdata.csv") #write out this output
  CombineReport[[1]] <- dat_rank_all
  
  datWide_Area <- dat_rank_all  %>%
    
    select(SampID,  MassFeature_Column, Adjusted_Area_VolNormed) %>%
    spread(., SampID, Adjusted_Area_VolNormed) %>%
    as.data.frame()
  ######   write_csv(datWide_Area, "WideArea.csv")
  CombineReport[[2]] <- datWide_Area
  
  
  datWide_Rank  <- dat_rank_all  %>%
    mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
    select(SampID, MassFeature_Column, Rank) %>%
    spread(., SampID, Rank) %>%
    as.data.frame()
##  write_csv(datWide_Rank, "WideRank.csv")
  CombineReport[[3]] <- datWide_Rank
  
  
  datWide_RankPercent  <- dat_rank_all  %>%
    mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
    select(SampID, MassFeature_Column, RankPercent) %>%
    spread(., SampID, RankPercent) %>%
    as.data.frame()
  CombineReport[[4]] <- datWide_RankPercent
 # write_csv(datWide_RankPercent, "WideRankPercentile.csv")
  
        
 return(CombineReport) 
}


ID_MFs <- function(id.file, dat.file, id.manual.file ){
  
  IDsetc <- read.csv(id.file) %>%
    select(MF_Frac, AvePoo:BestMatch)
  
  MFs <-  read.csv(dat.file) %>%
    select(MassFeature_Column) %>%
    separate(MassFeature_Column, into = c("MF_Frac", "Column"), sep = "_X_") %>%
    left_join(IDsetc)
  
  MFs2 <- MFs %>%
    select(MF_Frac:rt, BestMatch, Confidence, everything()) %>%
    mutate(BestMatch = ifelse(str_detect(MF_Frac, "^I\\d"), as.character(BestMatch), as.character(MF_Frac))) %>%
    mutate(Confidence = ifelse(str_detect(MF_Frac, "^I\\d"), Confidence, 1)) 
  
  #Edit them with the manual edits
  MFs_Fix <- read_csv(id.manual.file) 
  
  MFs_Fix2 <-  MFs_Fix %>%
    select(MF_Frac, Manual_Annote, Manual_Confidence) %>%
    filter(!is.na(Manual_Annote))
  
  MFs3 <- MFs2 %>%
    left_join(MFs_Fix2) %>%
    mutate(BestMatch = ifelse(is.na(Manual_Annote), as.character(BestMatch), as.character(Manual_Annote))) %>%
    mutate(Confidence = ifelse(is.na(Manual_Annote), Confidence, Manual_Confidence))
  
  datWide_RankPercent2 <- read.csv(dat.file) %>%
    separate(MassFeature_Column, into = c("MF_Frac", "Column"), sep = "_X_") %>%
    left_join(MFs3 %>% select(MF_Frac, BestMatch)) %>%
    mutate(MF_Frac = ifelse(is.na(BestMatch), MF_Frac, BestMatch)) %>%
    select(-BestMatch, -Column)
  

  return(datWide_RankPercent2)

}
  
# dat.file <-"Intermediates/WideArea_ModID.csv" 
# id.manual.file <- "RawOutput/MFs_match_wManual.csv"
# std.url = "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"

write_wide_w_MFinfo <- function(dat.file, id.manual.file, std.url ){

library(tidyverse)
library(RCurl)

dat <- read_csv(dat.file)
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

return(dat.final)
}




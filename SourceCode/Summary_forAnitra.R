#Make a table for Anitra of the culture data
cul.dat <- read_csv("Intermediates/Culture_Intermediates/combined_long_withquan.csv") %>%
  mutate(Org_Name = Org_Name %>% str_replace_all("CCMP22090", "") %>%
           str_replace_all("Emiliania huxleyi \nCCMP2090", "Emiliania huxleyi CCMP2090"))

mf.dat <- read_csv("Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv") %>%
  mutate(Cmp_ID = paste(Identification, Column, sep = "_"))%>%
  select(MassFeature_Column, Cmp_ID, Confidence)

cul.dat.mf.info <- cul.dat %>% left_join(mf.dat, by = "MassFeature_Column")%>%
  filter(Confidence == 1) %>%
  filter(!str_detect(ID_rep, "P3"))%>%
  filter(!str_detect(ID_rep, "P55"))

#Get biovolume normalized areas------
cul.dat.wide.full.bionorm <- cul.dat.mf.info %>% select(Cmp_ID, ID_rep, biovolArea) %>% unique() %>%
  spread(Cmp_ID, biovolArea)  %>%
  left_join(cul.dat %>% select(ID_rep, Org_Name, Org_Type) %>% unique(), by = "ID_rep") %>%
  select(Org_Type, ID_rep, Org_Name, everything())

cul.dat.wide.full.bionorm.mean <- cul.dat.mf.info %>% select(ID, Cmp_ID, biovolArea) %>%
  group_by(ID, Cmp_ID) %>%
  summarise(mean.biovolArea = mean(biovolArea, na.rm = TRUE)) %>%
  mutate(mean.biovolArea = ifelse(is.na(mean.biovolArea), NA, mean.biovolArea)) %>%
  spread(Cmp_ID, mean.biovolArea)  %>%
  left_join(cul.dat %>% select(ID, Org_Name, Org_Type) %>% unique(), by = "ID") %>%
  select(Org_Type, ID, Org_Name, everything()) 

#Get intracellular concentrations -----
cul.dat.wide.full.intracellconc <- cul.dat.mf.info %>% select(Cmp_ID, ID_rep, intracell_conc_umolCL) %>% unique() %>%
  spread(Cmp_ID, intracell_conc_umolCL)  %>%
  left_join(cul.dat %>% select(ID_rep, Org_Name, Org_Type) %>% unique(), by = "ID_rep") %>%
  select(Org_Type, ID_rep, Org_Name, everything())

cul.dat.wide.full.intracell.mean <- cul.dat.mf.info %>% select(ID, Cmp_ID, intracell_conc_umolCL) %>%
  group_by(ID, Cmp_ID) %>%
  summarise(mean.intracell = mean(intracell_conc_umolCL, na.rm = TRUE)) %>%
  mutate(mean.intracell = ifelse(is.na(mean.intracell), NA, mean.intracell)) %>%
  spread(Cmp_ID, mean.intracell)  %>%
  left_join(cul.dat %>% select(ID, Org_Name, Org_Type) %>% unique(), by = "ID") %>%
  select(Org_Type, ID, Org_Name, everything()) 


write_csv(cul.dat.wide.full.bionorm, "Intermediates/ForAnitra/Fulldat_AreaPerBiovol.csv")
write_csv(cul.dat.wide.full.bionorm.mean, "Intermediates/ForAnitra/Meandat_AreaPerBiovol.csv")
write_csv(cul.dat.wide.full.intracellconc, "Intermediates/ForAnitra/Fulldat_IntracellCConc.csv")
write_csv(cul.dat.wide.full.intracell.mean, "Intermediates/ForAnitra/Meandat_IntracellCConc.csv")

# id.file <- "RawOutput//NoIsos_allFractions_VolNormed_BMISd_CVFiltered_xsetWIDE_wMS2_wIDs.csv"
# id.manual.file <- "RawOutput/MFs_match_wManual.csv"
# std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"
# dat.hilicpos.file <- "RawOutput/HILICPos_IntegrationsBigPeaksWTargeted2.csv"
# dat.hilicneg.file <- "RawOutput/HILICNeg_IntegrationsBigPeaksWTargeted.csv"
# dat.cyano.file <- "RawOutput/CyanoAq_IntegrationsBigPeaksWTargeted.csv"
# samp.key.file <- "MetaData/Sample_key.csv"
# dat.file <- "Intermediates/Longdata.csv"
# is.names.file <- "MetaData/InternalStandardNames.csv"
# meta.data.file <- "MetaData/SampInfo_wMetaData.csv"
# 

#quan.dat.file <- "Intermediates/Quantified_LongDat.csv"
#metaG1.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"


###FUNCTION CALL GOES HERE
CMAPPER <- function(quan.dat.file, metaG1.file){

QuanDat <- read_csv(quan.dat.file)
MetaG1Dat <- read_csv(metaG1.file) %>%
    select(SampID, UTC_real, latitude, longitude, Depth) %>%
    mutate(Time = as.character(UTC_real)) %>% select(-UTC_real)
  
  quan.dat.withUTC <- left_join(QuanDat, MetaG1Dat) %>%
    mutate(pmolinEnviroave = nmolinEnviroave*1000)
  
  quan.dat.withUTC.wide.G1 <- quan.dat.withUTC %>%
    select(BestMatch, pmolinEnviroave, SampID, latitude, longitude, Time, Depth) %>%
    spread(key = BestMatch, value = pmolinEnviroave) %>%
    filter(!str_detect(SampID, "KM1513")) %>%
    filter(!str_detect(SampID, "S7E41"))
  
  quan.dat.withUTC.wide.G2 <- quan.dat.withUTC %>%
    select(BestMatch, pmolinEnviroave, SampID, latitude, longitude, Time, Depth) %>%
    spread(key = BestMatch, value = pmolinEnviroave) %>%
    filter(!str_detect(SampID, "KM1513")) %>%
    filter(str_detect(SampID, "S7E41"))
  
  CMAP <- list()
  CMAP[[1]] <- quan.dat.withUTC.wide.G1
  CMAP[[2]] <- quan.dat.withUTC.wide.G2
  
  return(CMAP)

}


#write.csv(quanDatSum, "Quantified_MFSummary.csv")
#write.csv(quanDat4, "Quantified_LongDat.csv")
#write.csv(quanDatWide, "Quantified_WideDat.csv")


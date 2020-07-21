#TO DO: take in the WideArea_withIDinfo_withCultureLogBioArea.csv file; get MS2s from the environment, write out a DF of the MS2 and what sample its from
library(tidyverse)
library(xcms)

#Name your inputs
dat.filename<- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
HILICPos.MS2.dir <- "../HILICPos_mzXML/DDA"
HILICNeg.MS2.dir <- "../HILICNeg_mzXML/DDA"
RP.MS2.dir <- "../CyanoAq_mzXML/DDA"


#Define the getms2spectra function---------
getms2spectra <-
  function(xs=xs, mass=mass, time=time, timetol=15, masstolLR=0.4, masstolHR=0.02){
    #xs is msn2xcmsRaw(xcmsRaw(DDAFILE, includeMSn=TRUE))
    masslow<-mass-masstolLR
    masshigh<-mass+masstolLR
    masslowhr <- mass - masstolHR
    masshighr <- mass + masstolHR
    timelow<-time-timetol
    timehigh<-time+timetol
    scanlist<-which(xs@msnPrecursorMz>masslow & xs@msnPrecursorMz<masshigh & xs@scantime>timelow & xs@scantime<timehigh)
    allMS2list <- list()
    if(length(scanlist)>0){
      for (l in 1:length(scanlist)){
        ms2dat <- as.data.frame(getScan(xs,scanlist[l])) %>% mutate(scannumber = scanlist[l])
        allMS2list[[l]] <- ms2dat} #This loop gets all matching scans with LR
      allMS2_matched <- do.call(rbind, allMS2list) %>%
        filter(mz > masslowhr & mz < masshighr) 
      matchedscanlist <- allMS2_matched$scannumber} else {matchedscanlist <- c()} #Do we find any HR matches?  
    if(length(matchedscanlist > 0)){
      MS2Summed <- do.call(rbind, allMS2list) %>%
        filter(scannumber %in% allMS2_matched$scannumber) %>%
        mutate(mzRND = as.factor(round(mz, digits = 2))) %>%
        group_by(mzRND) %>%
        summarise(mz = mean(mz),
                  intensity = sum(intensity),
                  sqrtintensity = sqrt(intensity),
                  scancount = n())
      bestscan <- MS2Summed %>%
        arrange(desc(intensity)) %>%
        mutate(intensity = round(intensity/max(intensity)*100, digits = 1))%>%
        filter(intensity > 0.5) %>%
        mutate(mz = round (mz, digits = 5),
               mash = paste(mz, intensity, sep = ", " ))
      sortedscanrange <- paste(bestscan$mash, collapse = "; ")
      return(sortedscanrange)}else{return(NA)} #Close if/then for scanlist >1
    
  }
#Close Function



#Get the data into the three separate fractions----
dat.all <- read_csv(dat.filename)
dat.HILICpos <- dat.all %>% filter(Column == "HILIC" & z == 1)
dat.HILICneg <- dat.all %>% filter(Column == "HILIC" & z == -1)
dat.RP <- dat.all %>% filter(Column == "RP")


#First go and get HILICPos MS2-----
dat.files.HILICpos <- list.files(HILICPos.MS2.dir, recursive = TRUE, full.names = TRUE) 
dat.files.HILICpos <- dat.files.HILICpos[grepl(".mzXML", dat.files.HILICpos)]

dat.HILICpos2 <- dat.HILICpos %>%
  select(MassFeature_Column, mz, rt)
names <- c("MS2dat1", "MS2dat2", "MS2dat3", "MS2dat4", "MS2dat5")

timetol = 20

for (k in 1:length(dat.files.HILICpos)){
  print(paste("Extracting MS2 from file", k, "of", length(dat.files.HILICpos), sep = " "))
  xs <- msn2xcmsRaw(xcmsRaw(dat.files.HILICpos[k], includeMSn=TRUE))
  dat.HILICpos2[, names[k]] <- NA %>% as.character()
  
  for (q in 1:length(dat.HILICpos2$mz))#Loop thru each MF
  {
    mass = dat.HILICpos2$mz[q]
    time = dat.HILICpos2$rt[q]
    ms2.spec<- getms2spectra(xs=xs, mass=mass, time=time, timetol=timetol)
    dat.HILICpos2[q, names[k]] <- ms2.spec
    
  } #Close MF loop
  
} #Close ms2 data loop

#Second go and get HILICNeg MS2-----
dat.files.HILICneg <- list.files(HILICNeg.MS2.dir, recursive = TRUE, full.names = TRUE) 
dat.files.HILICneg <- dat.files.HILICneg[grepl("\\d.mzXML", dat.files.HILICpos)]

dat.HILICneg2 <- dat.HILICneg %>%
  select(MassFeature_Column, mz, rt)
names <- c("MS2dat1", "MS2dat2", "MS2dat3", "MS2dat4", "MS2dat5")

timetol = 20

for (k in 1:length(dat.files.HILICneg)){
  print(paste("Extracting MS2 from file", k, "of", length(dat.files.HILICneg), sep = " "))
  xs <- msn2xcmsRaw(xcmsRaw(dat.files.HILICneg[k], includeMSn=TRUE))
  dat.HILICneg2[, names[k]] <- NA %>% as.character()
  
  for (q in 1:length(dat.HILICneg2$mz))#Loop thru each MF
  {
    mass = dat.HILICneg2$mz[q]
    time = dat.HILICneg2$rt[q]
    ms2.spec<- getms2spectra(xs=xs, mass=mass, time=time, timetol=timetol)
    dat.HILICneg2[q, names[k]] <- ms2.spec
    
  } #Close MF loop
  
} #Close ms2 data loop

#Third go and get RP MS2-----
dat.files.RP <- list.files(RP.MS2.dir, recursive = TRUE, full.names = TRUE) 
dat.files.RP <- dat.files.RP[grepl("\\d.mzXML", dat.files.HILICpos)]

dat.RP2 <- dat.RP %>%
  select(MassFeature_Column, mz, rt)
names <- c("MS2dat1", "MS2dat2", "MS2dat3", "MS2dat4", "MS2dat5")

timetol = 10

for (k in 1:length(dat.files.RP)){
  print(paste("Extracting MS2 from file", k, "of", length(dat.files.RP), sep = " "))
  xs <- msn2xcmsRaw(xcmsRaw(dat.files.RP[k], includeMSn=TRUE))
  dat.RP2[, names[k]] <- NA %>% as.character()
  
  for (q in 1:length(dat.RP2$mz))#Loop thru each MF
  {
    mass = dat.RP2$mz[q]
    time = dat.RP2$rt[q]
    ms2.spec<- getms2spectra(xs=xs, mass=mass, time=time, timetol=timetol)
    dat.RP2[q, names[k]] <- ms2.spec
    
  } #Close MF loop
  
} #Close ms2 data loop







#Now combine them back together
dat.ms2.all <- full_join(dat.HILICneg2, dat.HILICpos2) %>%
  full_join(dat.RP2)

dat.ms2.all.2 <- dat.ms2.all %>%
  pivot_longer(cols = MS2dat1:MS2dat4, 
               names_to = "MSdat", values_to = "Spectra") %>%
  filter(!is.na(Spectra)) %>%
  group_by(MassFeature_Column, mz, rt) %>%
  summarise(MS2 = Spectra[1])

dat.all.2 <- dat.all %>%
  select(MassFeature_Column:z) %>%
  left_join(dat.ms2.all.2, by = c("MassFeature_Column", "mz", "rt") )

write_csv(dat.all.2, "MF_info_withMS2s.csv")
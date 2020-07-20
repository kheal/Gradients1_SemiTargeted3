#TO DO: take in the WideArea_withIDinfo_withCultureLogBioArea.csv file; get MS2s from the environment, write out a DF of the MS2 and what sample its from
library(tidyverse)
library(xcms)

#Name your inputs
dat.filename<- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
HILICPos.MS2.dir <- "../HILICPos_mzXML/DDA"

#Get the data into the three fractions
dat.all <- read_csv(dat.filename)
dat.HILICpos <- dat.all %>% filter(Column == "HILIC" & z == 1)

#First go and get HILICPos 
dat.files.HILICpos <- list.files(HILICPos.MS2.dir, recursive = TRUE, full.names = TRUE) 
dat.files.HILICpos <- dat.files.HILICpos[grepl(".mzXML", dat.files.HILICpos)]

dat.HILICpos2 <- dat.HILICpos %>%
  select(MassFeature_Column, mz, rt)
names <- c("MS2dat1", "MS2dat2", "MS2dat3", "MS2dat4", "MS2dat5")

for (k in 1:length(dat.files.HILICpos)){
  print(paste("Extracting MS2 from file", k, "of", length(dat.files.HILICpos), sep = " "))
  xs <- msn2xcmsRaw(xcmsRaw(dat.files.HILICpos[k], includeMSn=TRUE))
  MZdat[, names[k]] <- NA 
  
  for (q in 1:length(MZdat$mz))#Loop thru each MF
  {
    mass = MZdat$mz[q]
    time = MZdat$rt[q]
    MZdat[q, names[k]] <- getms2spectra(xs=xs, mass=mass, time=time, timetol=timetol)
  } #Close MF loop
  
} #Close ms2 data loop

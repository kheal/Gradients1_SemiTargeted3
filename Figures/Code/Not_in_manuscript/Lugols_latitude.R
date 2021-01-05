library(tidyverse)
library(cowplot)
library(here)
library(plotly)
library(Cairo)
library(limSolve)

#Name your inputs
lugols.dat.file <- "MetaData/Lugols/Gradients_2016_Lugols Stations_format.csv"

#Plot PC and PN vs latitude - maybe add this to the map figure?
dat.all <- read_csv(lugols.dat.file) 

dat.long <- dat.all %>%
  gather(key = "org", value = "cellpermL", centric:protist)

#All separate
c <- ggplot(data = dat.long, aes(x =latitude, y =  cellpermL, color = org))+
  geom_point(size = 2)+
  facet_wrap(~ org, scales = "free")

ggsave("Figures/Preliminary/lugols.pdf", c, width = 10.5, height = 6, units = "in")

source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
library(here)
library(plotly)
library(beyonce)


#Name your inputs
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
long.dat.file <- "Intermediates/Longdata.csv"
cruise.id1 <- "KOK1606"


#Plot number 1 - All the enviro samples, colored by depth, shaped by cruise -----
#Read in long dat, calculated PC normed Adjusted_area, mudge to get a matrix of all the enviro samples
meta.dat.sub <- read_csv(here(meta.dat.file)) %>%
  select(SampID, PC_ave)

long.dat.all <- read_csv(long.dat.file)
long.dat.sub <- long.dat.all %>% select(MassFeature_Column, SampID, Adjusted_Area_VolNormed) %>%
  left_join(meta.dat.sub) %>%
  mutate(area.pc.normed = Adjusted_Area_VolNormed/PC_ave) %>%
  select(-Adjusted_Area_VolNormed, -PC_ave)

#Make data wide
combo.wide <- long.dat.sub %>%
  spread(key = SampID, value = area.pc.normed)
meta.dat <- read_csv(here(meta.dat.file)) %>% 
  filter(SampID %in% colnames(combo.wide))
enviro.names <- meta.dat$SampID

# Make the combo.file into a matrix
combo.wide.matrix<- combo.wide %>% select(c(enviro.names)) %>% as.data.frame()
row.names(combo.wide.matrix) <- combo.wide$MassFeature_Column


#Plot number 3 - Only KOK1601, colored by latitude or temperature
#Read in long dat, mudge to get a matrix of all the enviro samples
meta.dat.sub2 <- read_csv(here(meta.dat.file)) %>% 
  filter(SampID %in% colnames(combo.wide)) %>%
  filter(Cruise == "KOK1606")
enviro.names.sub2 <- meta.dat.sub2$SampID

# Make the combo.file into a matrix
combo.wide.matrix.sub2<- combo.wide %>% select(c(enviro.names.sub2)) %>% as.data.frame()
row.names(combo.wide.matrix.sub2) <- combo.wide$MassFeature_Column

# Standardize by z scores ### 
datwidestd.sub2 <- data.stand((combo.wide.matrix.sub2), method='standardize', margin='row', plot=F)

#WVdat.t.std <- WVdat.t.std[, !colSums(!is.finite(as.matrix(WVdat.t.std)))]
nmds.sub2<-metaMDS(t(datwidestd.sub2), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)

#vec.tp<-envfit(tp.nmds$points, TPExpDat.std, perm=1000)
pal <- rev(beyonce_palette(41, 160, type = "continuous"))

pointlocation.sub2 <- nmds.sub2[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.sub2[['points']])) %>%
  left_join(meta.dat) %>%
  mutate(lat_round = round(latitude, digits = 1))

c<- ggplot(data = pointlocation.sub2, aes(x =MDS1, y =  MDS2, group = lat_round,
                                          colour = latitude,
                                          shape = AMorPM,
                                          fill = latitude,
                                           label = SampID))+
  scale_fill_gradientn(colours = pal,
                       guide = guide_legend(nrow=10, keyheight = .5, keywidth = .5 ))+
  scale_colour_gradientn(colours = pal,
                       guide = guide_legend(nrow=10, keyheight = .5, keywidth = .5 ))+
  geom_polygon(fill = NA, color = "grey") +
  geom_point(size = 5) 
c  
d<- ggplot(data = pointlocation.sub2, aes(x =MDS1, y =  MDS2, group = lat_round,
                                          colour = temp1,
                                          shape = NBorSB,
                                          fill = temp1,
                                          label = SampID))+
  scale_fill_gradientn(colours = pal,
                       guide = guide_legend(nrow=10, keyheight = .5, keywidth = .5 ))+
  scale_colour_gradientn(colours = pal,
                         guide = guide_legend(nrow=10, keyheight = .5, keywidth = .5 ))+
  geom_polygon(fill = NA, color = "grey") +
  geom_point(size = 5) 
d  

save_plot("Figures/Manuscript_figures/NMDS_area_PCnormed.pdf", d, base_height = 6.5, base_width = 6.5)


  


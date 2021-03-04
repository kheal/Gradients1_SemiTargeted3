source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(beyonce)

#TO DO: save out the vectors with appropriate names?  Maybe just add these pvalues to the actual MF info data file

#Name inputs-----
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
combo.wide.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
cruise.id1 <- "KOK1606"

#Read in long dat, mudge to get a matrix of all the enviro samples----
long.dat.all <- read_csv(combo.wide.file)
combo.wide <- read_csv(combo.wide.file)
meta.dat <- read_csv(meta.dat.file) %>% 
  filter(SampID %in% colnames(combo.wide))
enviro.names <- meta.dat$SampID

# Make the combo.file into a matrix-----
combo.wide.matrix<- combo.wide %>% dplyr::select(c(enviro.names)) %>% as.data.frame()
row.names(combo.wide.matrix) <- combo.wide$MassFeature_Column

#Plot only KOK1606, colored by latitude or temperature-----
meta.dat.sub2 <- read_csv(meta.dat.file) %>% 
  filter(SampID %in% colnames(combo.wide)) %>%
  filter(Cruise == "KOK1606")
enviro.names.sub2 <- meta.dat.sub2$SampID

# Make the combo.file into a matrix, replace all NAs with 0s
combo.wide.matrix.sub2<- combo.wide %>% 
  dplyr::select(c(enviro.names.sub2)) %>% 
  as.data.frame()
row.names(combo.wide.matrix.sub2) <- combo.wide$MassFeature_Column
combo.wide.matrix.sub2[is.na(combo.wide.matrix.sub2)] <- 0

# Standardize by z scores, run NMDS, extract point location----
datwidestd.sub2 <- data.stand((combo.wide.matrix.sub2), method='standardize', margin='row', plot=F)
nmds.sub2<-metaMDS(t(datwidestd.sub2), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue <-nmds.monte(t(datwidestd.sub2), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result <- monte.pvalue[[2]]
print(paste(monte.pvalue.result, "= pvalue of nmds"))

pal <- rev((beyonce_palette(41, 160, type = "continuous")))
pointlocation.sub2 <- nmds.sub2[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.sub2[['points']])) %>%
  left_join(meta.dat, by = "SampID") %>%
  mutate(lat_round = round(latitude, digits = 1)) %>%
  mutate(Zone = factor(ifelse(Zone == 1, 1, 2)))

#Plot out the point location for the NMDS----
d<- ggplot(data = pointlocation.sub2, aes(x =MDS1, y =  MDS2, group = lat_round,
                                          colour = latitude,
                                          fill = latitude,
                                          label = SampID,
                                          shape = Zone))+
  scale_fill_gradientn(colours = pal, limits = c(23,38), breaks = c(24, 28, 32, 36), 
                       name = "Latitude")+
  scale_colour_gradientn(colours = pal, limits = c(23,38), breaks = c(24, 28, 32, 36), 
                         name = "Latitude")+
  annotate("text", x = -12, y = 12, 
           label = paste0("Stress = ", round(nmds.sub2[['stress']], digits = 2)), size = 2.2)+
  geom_polygon(fill = NA, color = "grey") +
  geom_point(size = 3) +
  guides(shape = FALSE)+
  theme(axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))

d  

#save_plot("Figures/Manuscript_figures/NMDS_area.pdf", d, base_height = 4.25, base_width = 4.5)


#Get and write out vector info, write out CSV into the correct folder Not sure this is necessary right now,  it's slow, but it works.  Commenting out for the time being.-----
#vec.info<-envfit(nmds.sub2$points, t(datwidestd.sub2), perm=9999)
#names <- as.data.frame(vec.info[[1]][[1]])%>%
#  mutate(Compound.Name = row.names(.),
#         pvalues = vec.info[[1]][[4]]) %>%
#  select(Compound.Name, pvalues)

#Assign zone info---- (NPSG = 1, NPTZ = 2)
samp.info <- as.data.frame(row.names(t(datwidestd.sub2)))
colnames(samp.info) <- "SampID"
samp.info <- samp.info %>%
  left_join(meta.dat.sub2 %>% dplyr::select(SampID, Zone, AMorPM, NBorSB, SampType), by = "SampID") %>%
  mutate(Zone = ifelse(Zone == 3, 2, Zone))

#Run ANOSIMS and save out 
Zone.anosim<-anosim(t(datwidestd.sub2),samp.info$Zone, distance = 'euclidean')
Zone.anosim.Results <- c("NPSG vs NPTZ", Zone.anosim$statistic, Zone.anosim$signif)

Day.anosim<-anosim(t(datwidestd.sub2),samp.info$AMorPM, distance = 'euclidean')
Day.anosim.Results <- c("Day vs Night", Day.anosim$statistic, Day.anosim$signif)

Direction.anosim<-anosim(t(datwidestd.sub2),samp.info$NBorSB, distance = 'euclidean')
Direction.anosim.Results <- c("North vs South transect", Direction.anosim$statistic, Direction.anosim$signif)

Samp.method.anosim<-anosim(t(datwidestd.sub2),samp.info$SampType, distance = 'euclidean')
Samp.method.anosim.Results <- c("Sampling method", Samp.method.anosim$statistic, Samp.method.anosim$signif)

AnosimResults <- rbind(Zone.anosim.Results, Day.anosim.Results) %>%
  rbind(Direction.anosim.Results) %>% as_tibble() 
colnames(AnosimResults) <- c("Variable Tested", "ANOSIM statistic", "p-value")
AnosimResults <- AnosimResults %>% arrange(`p-value`) %>%
  mutate(`ANOSIM statistic` = round(as.numeric(`ANOSIM statistic`), digits = 3))

write_csv(AnosimResults, "Intermediates/ANOSIM_results")



  


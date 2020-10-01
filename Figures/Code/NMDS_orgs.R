source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)

#Name inputs-----
dat.file <- "Intermediates/Culture_Intermediates/combined_long_cultures.csv"
meta.dat.file <- "MetaData/CultureMetaData.csv"

#Read in long dat, toss 32ppt samples, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
meta.dat <- read_csv(meta.dat.file) %>%
  rename(ID_rep = CultureID)

long.dat <- read_csv(dat.file) 
wide.dat <- long.dat %>%
  pivot_wider(id_cols = MassFeature_Column, names_from = ID_rep, values_from = LogBioArea)
wide.matrix<- wide.dat %>% dplyr::select(-MassFeature_Column) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$MassFeature_Column
compound.all.zeros <- wide.dat %>%
  dplyr::select(MassFeature_Column) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)
wide.matrix.2 <- wide.matrix[compound.all.zeros$MassFeature_Column, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0

#Run NMDS with only row standaridation, extract point location
#Stress is (nearly) zero: you may have insufficient data; samples are too different :(
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='max', margin='row', plot=F)
nmds.raw<-metaMDS(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(ID_rep = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "ID_rep") 

#Plot out the point location for the raw NMDS----
pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(8)[1:8])
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, 
                                              colour = Org_Type_Specific,
                                              fill = Org_Type_Specific, label = Org_Name))+
  geom_point(size = 3) +
  geom_polygon(aes(group = CultureID_short, color = Org_Type_Specific, fill= Org_Type_Specific), 
               alpha = 0.2) +
  scale_fill_manual(values = pal)+
  scale_color_manual(values = pal)+
  annotate("text", x = -3.4, y = 2.8, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 2), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 2.2)+
    theme(axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_blank(), 
        aspect.ratio=1)

d.raw  

print(d.raw)
save_plot("Figures/Manuscript_figures/NMDS_Organisms.pdf", d.raw, base_height = 5, base_width = 6)

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
PA.file <- "Intermediates/PA_table_forK.csv"

#Read in PA dat, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
PA.dat <- read_csv(PA.file) 
wide.matrix<- PA.dat %>% dplyr::select(-Organism, -Domain, -Group) %>% as.matrix()
row.names(wide.matrix) <- PA.dat$Organism
wide.matrix <- t(wide.matrix)
meta.dat <- PA.dat %>% dplyr::select(Organism, Domain, Group)

#Run NMDS, extract point location
nmds.raw<-metaMDS(t(wide.matrix), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(Organism = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "Organism") 


#Plot out the point location for the raw NMDS----
pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(8)[1:8])
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Group,
                                              colour = Group,
                                              fill = Group))+
  geom_point(position = position_jitter(width = 0.05, height = 0.05), size = 3) +
 scale_fill_manual(values = pal)+
  scale_color_manual(values = pal)+
  annotate("text", x = -3.4, y = 2.8, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 2), 
                          "\n p = ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 2.2)+
    theme(axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

d.raw  
print(d.raw)
save_plot("Figures/Preliminary/NMDS_PA.pdf", d.raw, base_height = 6, base_width = 6)

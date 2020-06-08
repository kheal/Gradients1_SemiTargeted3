library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(beyonce)

dat <- read_csv("Intermediates/Quantified_LongDat_Enviro.csv") %>% mutate(BestMatch = Identification)
dat2 <- dat %>%
  select(BestMatch, SampID,  nmolinEnviroave, nmolCave, percentCave, molFractionC, RankPercent )

datcul <- read_csv("Intermediates/culture_combined_long.csv") %>% 
  mutate(MassFeature_Column = MassFeature_Column %>% str_replace_all("HILICNeg", "HILIC") %>% str_replace_all("HILICPos", "HILIC")) %>%
           left_join(dat %>% select(BestMatch, MassFeature_Column) %>% unique)
dat2cul <- datcul %>%
  select(BestMatch, MassFeature_Column,ID_rep, PresAbs ) %>%
  filter(!is.na(BestMatch))

#Get a logical order for the MFs
AllSmps_Ordered <- dat %>%
  group_by(BestMatch) %>%
  summarise(AveSmp = median(nmolCave, na.rm = T)) %>%
  arrange((AveSmp)) %>%
  filter(!is.na(AveSmp))

#Get metadata to plot from both Skyline and XCMS
Metadatall <-  read_csv("MetaData/SampInfo_wMetaData.csv")
Metadat_rear <- Metadatall %>%
  arrange(latitude) 

#Get environmental data into a plottable shape
AllSmps_long <- dat2 %>%
  mutate(SampID = factor(SampID, levels= c(Metadat_rear$SampID))) %>%
  mutate(MF_Frac = factor(BestMatch, levels= c(AllSmps_Ordered$BestMatch))) %>%
  filter(MF_Frac %in% AllSmps_Ordered$BestMatch) %>%
  left_join(Metadat_rear %>% select(SampID, Depth)) 

AllSmps_long_surface <- AllSmps_long %>% filter(Depth < 30)

EnviroSumm <- AllSmps_long_surface %>%
  mutate(Observed = ifelse((is.na(nmolCave)), 0, 1)) %>%
  group_by(MF_Frac) %>%
  summarise(Obvs = sum(Observed)/41) %>%
  mutate(Type = "Surface")

EnviroSumm_deep <- AllSmps_long %>% filter(Depth > 30) %>%
  mutate(Observed = ifelse((is.na(nmolCave)), 0, 1)) %>%
  group_by(MF_Frac) %>%
  summarise(Obvs = sum(Observed)/18) %>%
  mutate(Type = "Subsurface")

#Get culture data into a plottable shape
CultureSumm <- dat2cul %>%
  mutate(MF_Frac = factor(BestMatch, levels= c(AllSmps_Ordered$BestMatch))) %>%
  filter(MF_Frac %in% AllSmps_Ordered$BestMatch) %>%
  group_by(MF_Frac) %>%
  summarise(Obvs = sum(PresAbs)/66) %>%
  mutate(Type = "Culture")

CombinedSumm <- bind_rows(EnviroSumm, CultureSumm) %>% bind_rows(EnviroSumm_deep)

#Make a box plot of the range of concentrations in surface particles
pal <- rev(beyonce_palette(41, 11, type = "continuous"))

c <- ggplot(data = AllSmps_long, aes(factor(MF_Frac), nmolinEnviroave)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.3) +
  scale_y_log10()+  
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text (size = 7),
        axis.title.x= element_text (size = 8),
        panel.grid.major.y = element_line(colour="lightgrey", size = rel(0.2)))+
  labs(ylab("nomlC of quantified compounds"))

d <- ggplot(data = CombinedSumm, aes(factor(MF_Frac), y = Type, fill = Obvs)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal,
                      guide = guide_legend(reverse = TRUE, nrow=5, keyheight = .5, keywidth = .5 ),
                       limits = c(0,1))+
  coord_flip() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=-60, hjust=0, size = 7),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_blank())
 
f <- plot_grid(c, d,
                   align = "h", axis = "bt",
                   rel_widths = c(8, 2),
                   ncol = 2, scale = 1)



save_plot("Figures/Manuscript_figures/summaryBoxPlot_nmolCompound.pdf", f, base_height = 6.5, base_width = 6)






library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RColorBrewer)
library(viridis)
library(RCurl)

#Name your files -----
quandat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
culture.dat.long.filename <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
culture.meta.dat.filename <- "MetaData/CultureMetaData.csv"

#Get list of better names
std.url <- 'https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv'
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound_Name_Original,
         BestMatch =Compound_Name_Figure) %>%
  select(BestMatch, Identification) %>% unique()

#Load up the quan enviro data-----
dat <- read_csv(quandat.file) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)
dat2 <- dat %>%
  select(Identification, SampID,  nmolCave, molFractionC, RankPercent )

#Load up the culture data ----
dat2cul <- read_csv(culture.dat.long.filename) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)%>%
  select(Identification, MassFeature_Column,ID_rep, PresAbs, intracell_conc_umolCL) %>%
  filter(!is.na(Identification))

cult.meta.dat <- read_csv(culture.meta.dat.filename) %>% rename(ID_rep = CultureID) %>% select(ID_rep, Org_Type)
dat2cul <- dat2cul %>% left_join(cult.meta.dat, by = "ID_rep") %>%
  filter(!str_detect(ID_rep, "P3")) %>%
  filter(!str_detect(ID_rep, "P55"))


#Get a logical order for the MFs
AllSmps_Ordered <- dat %>%
  group_by(Identification) %>%
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
  mutate(MF_Frac = factor(Identification, levels= c(AllSmps_Ordered$Identification))) %>%
  filter(MF_Frac %in% AllSmps_Ordered$Identification) %>%
  left_join(Metadat_rear %>% select(SampID, Depth), by = "SampID") 

AllSmps_long_surface <- AllSmps_long %>% filter(Depth < 30)

dat3cul <- dat2cul%>% 
  mutate(Identification = factor(Identification, levels= c(AllSmps_Ordered$Identification))) %>%
  filter(!is.na(Identification))
  


#Summarize pres/abs in environmental samples and in culture samples
EnviroSumm <- AllSmps_long_surface %>%
  mutate(Observed = ifelse((is.na(nmolCave)), 0, 1)) %>%
  group_by(MF_Frac) %>%
  summarise(Obvs = sum(Observed)/n()) %>%
  mutate(Type = "Environmental") %>%
  mutate(Type_specific = "Surface")

EnviroSumm_deep <- AllSmps_long %>% filter(Depth > 30) %>%
  mutate(Observed = ifelse((is.na(nmolCave)), 0, 1)) %>%
  group_by(MF_Frac) %>%
  summarise(Obvs = sum(Observed)/n()) %>%
  mutate(Type = "Environmental") %>%
  mutate(Type_specific = "Subsurface")

#Get culture data into a plottable shape
CultureSumm <- dat2cul %>%
  mutate(MF_Frac = factor(Identification, levels= c(AllSmps_Ordered$Identification))) %>%
  filter(MF_Frac %in% AllSmps_Ordered$Identification) %>%
  group_by(MF_Frac, Org_Type) %>%
  summarise(Obvs = sum(PresAbs)/n()) %>%
  mutate(Type = "Culture") %>%
  rename(Type_specific = Org_Type)

CombinedSumm <- bind_rows(EnviroSumm, CultureSumm) %>% bind_rows(EnviroSumm_deep) %>%
  mutate(Type = factor(Type, levels = c("Environmental", "Culture"))) %>%
  mutate(Type_specific = factor(Type_specific, levels = c("Surface", "Subsurface" , "Cyanobacteria","Diatom", "Dinoflagellate", "Haptophyte","Prasinophyte", "Archaea"  )))

#Make a box plot of the range of concentrations in surface particles
pal <- rev(beyonce_palette(41, 11, type = "continuous"))

c <- ggplot(data = AllSmps_long, aes(factor(MF_Frac), nmolCave)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.3) +
  scale_y_log10(breaks = c(0.00001,  0.001,  0.1, 10, 1000), labels = function(x) sprintf("%g", x))+  
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text (size = 5),
        axis.title.x= element_text (size = 6, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        panel.grid.major.y = element_line(colour="lightgrey", size = rel(0.2)),
        plot.margin = unit(c(0,0,0,0), "mm"))+
  labs(y = "nmol C / L \nenvironmental concentrations")


#Make a box plot of the range of intracellular concentrations in cultures
c2 <- ggplot(data = dat3cul, aes(factor(Identification), intracell_conc_umolCL/1000)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.3) +
  scale_y_log10(breaks = c(0.00001,  0.001,  0.1, 10, 1000), labels = function(x) sprintf("%g", x))+  
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text (size = 5),
        axis.title.x= element_text (size = 6, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        panel.grid.major.y = element_line(colour="lightgrey", size = rel(0.2)),
        plot.margin = unit(c(0,0,0,0), "mm"))+
  labs(y = "mmol C / L \nintracellcular concentrations")



d <- ggplot(data = CombinedSumm, aes(x = Type_specific, y = factor(MF_Frac), fill = Obvs)) +
  geom_tile() +
  facet_grid(~ Type, scales = "free_x", space = "free_x")+
  scale_fill_viridis(option = 1, na.value="grey", 
                       breaks = c(0, 0.5, 1.0)) +
  guides(fill = guide_colourbar(barwidth = 3, barheight = 0.7, ticks.colour = "black"))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=-60, hjust=0, size = 5),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 5, margin = margin(t = -3, r = 0, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification="left",
        legend.margin=margin(-15,0,0,0),#T, r, b, l
        #legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0,4,0,0), "mm"),
        panel.spacing = unit(1, "mm")) 
d

f <- plot_grid(c, c2, d,
                   align = "h", axis = "bt",
                   rel_widths = c(10, 6, 3),
                   ncol = 3, scale = 1)



f

save_plot("Figures/Manuscript_figures/summaryBoxPlot.pdf", f, base_height = 7.5, base_width = 6)






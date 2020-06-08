library(tidyverse)
library(here)
library(cowplot)
require(RColorBrewer)
library(colorRamps)

num_toplot <- 19
dat <- read_csv("Intermediates/Quantified_LongDat.csv") %>% mutate(BestMatch = Identification)
dat2 <- dat %>%
  select(BestMatch, SampID, percentCave, nmolCave)

dat3 <- dat2 %>%
  group_by(SampID) %>%
  summarise(TotalpercentC_measured = sum(percentCave, na.rm = T)) %>% ungroup()

dat2 <- dat2 %>%
  left_join(dat3)

#Get a logical order for the MFs, pull out only the top 20, put in the colors  here?
AllSmps_Ordered <- dat %>%
  group_by(BestMatch) %>%
  summarise(AveSmp = mean(nmolCave, na.rm = T)) %>%
  arrange(desc(AveSmp)) %>%
  head(num_toplot) %>%
  mutate(Colors = colorRampPalette(brewer.pal(9,"Set1"))(num_toplot+1)[1:num_toplot])


#Get metadata to plot
Metadatall <-  read_csv("MetaData/SampInfo_wMetaData.csv")
Metadat_rear <- Metadatall %>%
  arrange(latitude)

#Get data into a plottable shape
AllSmps_long <- dat2 %>%
  mutate(SampID = factor(SampID, levels= c(Metadat_rear$SampID))) %>%
  mutate(MF_Frac = as.character(BestMatch))%>%
  filter(MF_Frac %in% AllSmps_Ordered$BestMatch)

AllSmps_long_Summary <- AllSmps_long %>%
  left_join(Metadatall %>% select(SampID, Station, Zone, Depth, Cruise, latitude)) %>%
  mutate(Latitude = as.factor(round(latitude, digits =1))) %>%
  group_by(Latitude, MF_Frac, Zone, Depth, Cruise) %>%
  summarise(percentCave = mean(as.numeric(percentCave), na.rm = T),
            totalAveSmp = mean(TotalpercentC_measured, na.rm = T))%>%
  ungroup() %>%
  mutate(Latitude = as.factor(Latitude))

Others <- AllSmps_long_Summary %>%
  group_by(Latitude, Zone, Depth, Cruise, totalAveSmp) %>% 
  summarise(measured = sum(percentCave, na.rm = T)) %>% ungroup() %>%
  mutate(percentCave = totalAveSmp-measured,
         MF_Frac = "all others") %>% select(-measured)

AllSmps_long_Summary_withOthers <- bind_rows(AllSmps_long_Summary, Others) %>%
  left_join(AllSmps_Ordered %>% rename(MF_Frac = BestMatch) %>% select(MF_Frac,  Colors)) %>%
  mutate(MF_Frac = factor(MF_Frac, levels= c(AllSmps_Ordered$BestMatch, "all others"))) %>%
  mutate(Colors = ifelse(is.na(Colors), "#999999", Colors)) %>%
  arrange(MF_Frac)

#Time to plot it up Transect ------
#colourCount <- length(unique(AllSmps_long_Summary_withOthers$MF_Frac)) # number of levels
AllSmps_long_Summary_Tran <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise == "KOK1606")

b <- ggplot(AllSmps_long_Summary_Tran, aes(x = Latitude, y = percentCave, fill = MF_Frac))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,8.5))+
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_Tran$Colors)),
                    guide = guide_legend(nrow=5, keyheight = 1, keywidth = 1 ))+
  labs(y= "percent C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position="bottom",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))+
  coord_flip() 

#Plotting up the North Depth Profile -----
AllSmps_long_Summary_nDP <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise ==  "MGL1704") %>%
  unique()

c <- ggplot(AllSmps_long_Summary_nDP, aes(x = as.numeric(Depth), y = percentCave, fill = MF_Frac))+
  scale_x_reverse() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,900))+
  geom_bar(stat = "identity", color = "black", width = 15, size = 0.2) +
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_nDP$Colors)),
                    guide = guide_legend(nrow=5,keyheight = 0.8, keywidth = 0.8 ))+
  labs(y= "percent C", x = "Depth")+
  guides(fill=FALSE) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))+
  coord_flip() 

#Plotting up the South Depth Profile -----
AllSmps_long_Summary_sDP <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise ==  "KM1513") %>%
  unique()

d <- ggplot(AllSmps_long_Summary_sDP, aes(x = as.numeric(Depth), y = percentCave, fill = MF_Frac))+
  scale_x_reverse() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,5))+
  geom_bar(stat = "identity", color = "black", width = 15, size = 0.2) +
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_Tran$Colors)),
                    guide = guide_legend(nrow=5,keyheight = 0.8, keywidth = 0.8 ))+
  labs(y= "percent C", x = "Depth")+
  guides(fill=FALSE) +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))+
  coord_flip() 

#Combine the plots 
c_d <- plot_grid(c, d, ncol = 2, labels = c('A', 'B'))
b_c_d <- plot_grid(c_d, b, ncol = 1, rel_heights = c(0.6, 1), labels = c('', 'C'))

#save_plot("Figures/Preliminary/nMoleC.pdf", b_c_d, base_height = 8, base_width = 6.5, units = "in")

#Try different orientation
c_d2 <- plot_grid(c, d, ncol = 1, labels = c('B', 'C'))
b_c_d2 <- plot_grid(b, c_d2, ncol = 2, labels = c('A', ''))
save_plot("Figures/Manuscript_figures/barplot_nmolC_PCnormed.pdf", b_c_d2, base_height = 6, base_width = 6.5, units = "in")



library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
require(RColorBrewer)
library(colorRamps)
library(RCurl)

#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         BestMatch = Compound.Name_figure) %>%
  select(BestMatch, Identification) %>% unique()

num_toplot <- 19
dat <- read_csv("Intermediates/Quantified_LongDat_Enviro.csv") %>% 
  left_join(stds.dat, by = "Identification") %>% 
  select(-Identification)
dat2 <- dat %>%
  select(BestMatch, SampID, nmolCave, totalCmeasured_nM) %>%
  filter(BestMatch != "DMSP")
 
#Get a logical order for the MFs, pull out only the top 20, put in the colors  here?
AllSmps_Ordered <- dat %>%
  group_by(BestMatch) %>%
  summarise(AveSmp = mean(nmolCave, na.rm = T),
            totalAveSmp = mean(totalCmeasured_nM, na.rm = T)) %>%
  arrange(desc(AveSmp)) %>%
  head(num_toplot) %>%
  mutate(Colors = colorRampPalette(brewer.pal(9,"Set1"))(num_toplot+1)[1:num_toplot])


#Get metadata to plot
Metadatall <-  read_csv("MetaData/SampInfo_wMetaData_withUTC.csv")
Metadat_rear <- Metadatall %>%
  arrange(latitude)

#Get data into a plottable shape
AllSmps_long <- dat2 %>%
  mutate(SampID = factor(SampID, levels= c(Metadat_rear$SampID))) %>%
  mutate(MF_Frac = as.character(BestMatch))%>%
  filter(MF_Frac %in% AllSmps_Ordered$BestMatch)

AllSmps_long_Summary <- AllSmps_long %>%
  left_join(Metadatall %>% select(SampID, Station, Zone, Depth, Cruise, latitude, Station_1)) %>%
  group_by(Station_1, MF_Frac, Zone, Depth, Cruise) %>%
  summarise(Latitude = round(mean(latitude), digits = 1),
            nmolCave = mean(as.numeric(nmolCave), na.rm = T),
            totalAveSmp = mean(totalCmeasured_nM, na.rm = T))%>%
  ungroup() %>%
  mutate(Latitude = as.factor(Latitude))

Others <- AllSmps_long_Summary %>%
  group_by(Latitude, Zone, Depth, Cruise, totalAveSmp) %>% 
  summarise(measured = sum(nmolCave, na.rm = T)) %>% ungroup() %>%
  mutate(nmolCave = totalAveSmp-measured,
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

b <- ggplot(AllSmps_long_Summary_Tran, aes(x = Latitude, y = nmolCave, fill = MF_Frac))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  annotate("segment", x = 4.5, xend = 4.5,y = 0, yend = 250, size = 1.5, linetype=2)+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0,250))+
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_Tran$Colors)),
                    guide = guide_legend(nrow=5, keyheight = 1, keywidth = 1 ,
                                         label.hjust = 0))+
  labs(y= "nMole C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6, margin = margin(l = -5)),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))+
  coord_flip() 


#Plotting up the North Depth Profile -----
AllSmps_long_Summary_nDP <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise ==  "MGL1704") %>%
  unique()

c <- ggplot()+
  geom_bar(data= AllSmps_long_Summary_nDP %>%  filter(Depth == "40"), 
           aes(x = as.numeric(Depth), y = nmolCave, fill = MF_Frac),
           stat = "identity", color = "black", width = 15, size = 0.2)+
  geom_bar(data= AllSmps_long_Summary_nDP %>%  filter(Depth != "40"), 
           aes(x = as.numeric(Depth), y = nmolCave, fill = MF_Frac),
           stat = "identity", color = "black", width = 15, size = 0.2)+
  scale_x_reverse() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1000))+
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_nDP$Colors)),
                    guide = guide_legend(nrow=5,keyheight = 0.8, keywidth = 0.8 ))+
  labs(y= "nMole C", x = "Depth (m)")+
  guides(fill=FALSE) +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))+
  coord_flip() 

#Plotting up the South Depth Profile -----
AllSmps_long_Summary_sDP <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise ==  "KM1513") %>%
  unique()

d <- ggplot(AllSmps_long_Summary_sDP, aes(x = as.numeric(Depth), y = nmolCave, fill = MF_Frac))+
  scale_x_reverse() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  geom_bar(stat = "identity", color = "black", width = 15, size = 0.2) +
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_Tran$Colors)),
                    guide = guide_legend(nrow=5,keyheight = 0.8, keywidth = 0.8 ))+
  labs(y= "nMole C", x = "Depth (m)")+
  guides(fill=FALSE) +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))+
  coord_flip() 

#Combine the plots 
c_d2 <- plot_grid(c, d, ncol = 1, labels = c('B', 'C'))
b_c_d2 <- plot_grid(b, c_d2, ncol = 2, labels = c('A', ''))

save_plot("Figures/Manuscript_figures/barplot_nmolC.pdf", b_c_d2, base_height = 6, base_width = 6.5, units = "in")  



##Make the same plots but with better set of font sizes....-----
#Time to plot it up Transect ------
#colourCount <- length(unique(AllSmps_long_Summary_withOthers$MF_Frac)) # number of levels
AllSmps_long_Summary_Tran <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise == "KOK1606")

b.big <- ggplot(AllSmps_long_Summary_Tran, aes(x = Latitude, y = nmolCave, fill = MF_Frac))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,250))+
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_Tran$Colors)),
                    guide = guide_legend(nrow=5, keyheight = 1, keywidth = 1 ))+
  labs(y= "nMole C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))+
  coord_flip() 
b.big

#Plotting up the North Depth Profile -----
AllSmps_long_Summary_nDP <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise ==  "MGL1704") %>%
  unique()

c.big <- ggplot(AllSmps_long_Summary_nDP, aes(x = as.numeric(Depth), y = nmolCave, fill = MF_Frac))+
  scale_x_reverse() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1000))+
  geom_bar(stat = "identity", color = "black", width = 15, size = 0.2) +
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_nDP$Colors)),
                    guide = guide_legend(nrow=5,keyheight = 0.8, keywidth = 0.8 ))+
  labs(y= "nMole C", x = "Depth (m)")+
  guides(fill=FALSE) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 9))+
  coord_flip() 

#Plotting up the South Depth Profile -----
AllSmps_long_Summary_sDP <- AllSmps_long_Summary_withOthers %>%
  filter(Cruise ==  "KM1513") %>%
  unique()

d.big <- ggplot(AllSmps_long_Summary_sDP, aes(x = as.numeric(Depth), y = nmolCave, fill = MF_Frac))+
  scale_x_reverse() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,150))+
  geom_bar(stat = "identity", color = "black", width = 15, size = 0.2) +
  scale_fill_manual(values = as.character(unique(AllSmps_long_Summary_Tran$Colors)),
                    guide = guide_legend(nrow=5,keyheight = 0.8, keywidth = 0.8 ))+
  labs(y= "nMole C", x = "Depth (m)")+
  guides(fill=FALSE) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))+
  coord_flip() 

#Combine the plots and save  out -----
c_d2.big <- plot_grid(c.big, d.big, ncol = 1, labels = c('', ''))
b_c_d2.big <- plot_grid(b.big, c_d2.big, ncol = 2, labels = c('', ''))
b_c_d2.big

cairo_pdf("Figures/Presentation_figures/barplot_nmolC.pdf", family="Arial Unicode MS", 6.5,6)
b_c_d2.big
dev.off()



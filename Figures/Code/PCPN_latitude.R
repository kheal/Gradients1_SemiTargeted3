library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)

#Name your inputs
pc.dat.file <- "MetaData/PCPN/KOK1606_PCPN_UW_Preliminary_OSU_KRH.csv"
quan.dat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"

#Plot PC and PN vs latitude - maybe add this to the map figure?
dat <- read_csv(pc.dat.file) %>%
  mutate(NBorSB = ifelse(row_number() < 132, "NB", "SB"))

dat.quan <- read_csv(quan.dat.file) %>%
  select(SampID, totalCmeasured_nM, totalNmeasured_nM) %>% unique() %>%
  left_join(read_csv(meta.dat.file) %>% select(latitude, SampID, Cruise, NBorSB, AMorPM)) %>%
  filter(Cruise == "KOK1606")

b.metabs <- ggplot()+
 # geom_point(data = dat, aes(x =Latitude, y =  PC), size = 1)+
  geom_bar(data = dat.quan, aes(x = latitude, y = totalCmeasured_nM), stat = "identity", width = 0.2, 
           position = position_dodge2(preserve = "single"))+
  scale_fill_manual(labels = c("Day", "Night"), values = c("grey","black")) +
  scale_y_continuous(bquote("Total quantifiable \n metabolites"), 
                     #limits = c(0, 10), 
                     expand = c(0,0))+
  scale_x_continuous(limits = c(23, 37.5))+
  theme(legend.position = c(0.05, 0.8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.box.background = element_rect(colour = "black"),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
b.metabs

c<- ggplot()+
  geom_point(data = dat, aes(x =Latitude, y =  PN), size = 1)+
  geom_bar(data = dat.quan, aes(x = latitude, y = totalNmeasured_nM/1000*10, fill = AMorPM), stat = "identity", width = 0.2, 
           position = position_dodge2(preserve = "single"))+
  scale_fill_manual(labels = c("Day", "Night"), values = c("grey","black")) +
  scale_y_continuous(bquote("Total particulate nitrogen (\U003BCM N L" ^-1*')'), limits = c(0, 2), expand = c(0,0),
                     sec.axis = sec_axis(~ (./10), name = bquote("Total quantifiable metabolites (\U003BCM N L" ^-1*')')))+
  scale_x_continuous(limits = c(23, 37.5))+
  labs(fill = "Time of collection for \n metabolite samples")+
  theme(legend.position = c(0.05, 0.8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.box.background = element_rect(colour = "black"),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
c

#save_plot("Figures/Manuscript_figures/PC_lat.pdf", b, base_height = 4, base_width = 6.5, device=cairo_pdf)
#save_plot("Figures/Manuscript_figures/PN_lat.pdf", c, base_height = 4, base_width = 6.5, device=cairo_pdf)


  


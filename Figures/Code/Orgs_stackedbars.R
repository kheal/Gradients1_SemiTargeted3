library(tidyverse)
library(here)

dat.filename <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
Meta.dat.file <- "MetaData/CultureMetaData.csv"

#TO DO: Lots of work to do in here...

#Get quan data in shape to plot
dat <- read_csv(dat.filename)
meta.dat <- read_csv(Meta.dat.file)

dat.mean <- dat %>% ungroup () %>%
  group_by(Identification, ID, Org_Name, Org_Type_Specific, Org_Type) %>%
  mutate(intracell_conc_umolCL = ifelse(is.na(intracell_conc_umolCL), 0 , intracell_conc_umolCL)) %>%
  summarise(intracell_conc_umolCL = mean(intracell_conc_umolCL))

dat.total <- dat.mean %>%
  group_by(ID) %>%
  summarise(Total_nmolCMeasured = sum(intracell_conc_umolCL))

dat.mean.molefraction <- dat.mean %>% ungroup () %>%
  left_join(dat.total, by = "ID") %>%
  mutate(molfraction = intracell_conc_umolCL/Total_nmolCMeasured)

#Get good organism order
order.of.org.df <- dat %>% ungroup %>% 
  select(Org_Type, Org_Type_Specific, ID) %>%
  unique() %>% arrange(Org_Type, Org_Type_Specific)

order.of.org <- order.of.org.df$ID

dat.mean.molefraction$ID <- factor(dat.mean.molefraction$ID, 
                                          levels = (order.of.org))

#Get good compound order
order.of.compounds <- dat.mean %>% ungroup %>% 
  group_by(Identification) %>%
  summarise(intracell_conc_umolCL = mean(intracell_conc_umolCL, na.rm = TRUE)) %>%
  arrange(desc(intracell_conc_umolCL))

dat.mean.molefraction$Identification <- factor(dat.mean.molefraction$Identification, 
                                   levels = (order.of.compounds$Identification))

compounds.to.highlight <- head(order.of.compounds$Identification, 10)

dat.mean.molefraction.withcolors <- dat.mean.molefraction %>%
  mutate(color = ifelse(Identification %in% compounds.to.highlight, Identification, "Other"))

#Plot to highlight top 19 in the samples?
#pal <- c("#9B445D", rep("grey", 84))
pal <- c(colorRampPalette(brewer.pal(9,"Dark2"))(20)[1:19], rep("grey", 84))

b.all <- ggplot(dat.mean.molefraction.withcolors, aes(x = ID, y = molfraction, fill = Identification))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  # geom_bar(stat = "identity", color = "black", size = 0.2, aes(fill = color))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
b.all
ggsave("Figures/Preliminary/stackedbar_org.pdf", b.all, width = 8, height = 6, units = "in")




#Plot to highlight glutamic acid
pal <- c("#9B445D", rep("grey", 84))

b.glu <- ggplot(dat.mean.molefraction.withcolors, aes(x = ID, y = molfraction, fill = Identification))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  # geom_bar(stat = "identity", color = "black", size = 0.2, aes(fill = color))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
b.glu

#ggsave("Figures/Presentation_figures/stackedbar_org_GLU.pdf", b.glu, width = 8, height = 6, units = "in")


#Plot to highlight gonyol
pal <- c(rep("grey", 2), "#BD6066", rep("grey", 81))

b.gon <- ggplot(dat.mean.molefraction.withcolors, aes(x = ID, y = molfraction, fill = Identification))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  # geom_bar(stat = "identity", color = "black", size = 0.2, aes(fill = color))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
b.gon

#Plot to highlight homarine
pal <- c(rep("grey", 3), "red", rep("grey", 82))

b.hom <- ggplot(dat.mean.molefraction.withcolors, aes(x = ID, y = molfraction, fill = Identification))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
 # geom_bar(stat = "identity", color = "black", size = 0.2, aes(fill = color))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
b.hom

#ggsave("Figures/Presentation_figures/stackedbar_org_HOMARINE.pdf", b.hom, width = 8, height = 6, units = "in")


#Plot to highlight Glucosylglycerol
pal <- c(rep("grey", 1), "#BD6066", rep("grey", 83))

b.glugly <- ggplot(dat.mean.molefraction.withcolors, aes(x = ID, y = molfraction, fill = Identification))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  # geom_bar(stat = "identity", color = "black", size = 0.2, aes(fill = color))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
b.glugly

#ggsave("Figures/Presentation_figures/stackedbar_org_Glucosylglycerol.pdf", b.glugly, width = 8, height = 6, units = "in")

#Plot to highlight DHPS
pal <- c(rep("grey", 17), "#3C8A9B", rep("grey", 83))

b.DHPS <- ggplot(dat.mean.molefraction.withcolors, aes(x = ID, y = molfraction, fill = Identification))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  # geom_bar(stat = "identity", color = "black", size = 0.2, aes(fill = color))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
b.DHPS

#ggsave("Figures/Presentation_figures/stackedbar_org_DHPS.pdf", b.DHPS, width = 8, height = 6, units = "in")





#Plot to highlight homarine
pal <- c(rep("grey", 3), "red", rep("grey", 82))

b.hom <- ggplot(dat.mean.molefraction.withcolors, aes(x = ID, y = molfraction, fill = Identification))+
  geom_bar(stat = "identity", color = "black", size = 0.2)+
  # geom_bar(stat = "identity", color = "black", size = 0.2, aes(fill = color))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))
b.hom






dat.homarine <- dat %>%
  filter(Identification == "Homarine") %>%
  left_join(meta.dat, by = "SampID")
  
dat.homarine.med <- dat.homarine %>%
  group_by(Identification, Station_1, Cruise, Depth) %>%
  summarise(nmolinEnviroave.med = mean(nmolinEnviroave),
            nmolinEnviroave.sd = sd(nmolinEnviroave),
            latitude = mean(latitude)) %>%
  mutate(nmolinEnviroave.sd = ifelse(is.na(nmolinEnviroave.sd), .2*nmolinEnviroave.med, nmolinEnviroave.sd))

g.homarine.KOK <- ggplot() +
  geom_ribbon(data = dat.homarine.med %>% filter(Cruise == "KOK1606"), 
              aes(x = latitude, ymin = nmolinEnviroave.med-nmolinEnviroave.sd,
                  ymax = nmolinEnviroave.med+nmolinEnviroave.sd), alpha = 0.4)+
    geom_point(data = dat.homarine %>% filter(Cruise == "KOK1606"), 
               aes(x = latitude, y = nmolinEnviroave)) +
  xlab("Latitude") +
  ylab("nmol / L")

g.homarine.KOK

#ggsave("Figures/Presentation_figures/homarineKOK.pdf", width = 8, height = 5, units = "in")
  
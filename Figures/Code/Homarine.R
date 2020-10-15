library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())
library(RCurl)
#library(magick)


#Name your compounds of interest
cmpds <- c("Homarine", "Trigonelline")
transform_factor = 50

#Name your files -----
quandat.file <- "Intermediates/Quantified_LongDat_Enviro.csv"
culture.dat.long.filename <- "Intermediates/Culture_Intermediates/combined_long_withquan.csv"
culture.meta.dat.filename <- "MetaData/CultureMetaData.csv"
field.meta.dat.filename <- "MetaData/SampInfo_wMetaData_withUTC.csv"
chromat.filename <- "RawOutput/Homarine_EICs.csv"

#Get list of better names-----
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         BestMatch = Compound.Name_figure) %>%
  select(BestMatch, Identification) %>% unique()

#Load up metadata-----
meta.dat.enviro <- read_csv(field.meta.dat.filename) %>%
  select(SampID, Volume, PC_ave, Station_1, Cruise, Depth, latitude)

meta.dat.culture <- read_csv(culture.meta.dat.filename)%>%
  select(CultureID, CultureID_short, BioVol_perFilter_uL, nmolC_filtered_final, Org_Type, Species)

#Load up chromat data -----
chromat.dat <- read_csv(chromat.filename, skip = 1)

#Load up culture data ----
cult.dat <- read_csv(culture.dat.long.filename)
cul.dat2 <- cult.dat %>%
  filter(Identification %in% cmpds) %>%
  select(Identification, CultureID_short, Org_Name:Org_Type_Specific, intracell_conc_umolL) %>%
  mutate(intracell_conc_umolL = ifelse(is.na(intracell_conc_umolL), 0, intracell_conc_umolL)) %>%
  group_by(Identification, CultureID_short, Org_Name, Org_Type, Org_Type_Specific) %>%
  summarise(intracell_conc_mmolL = mean(intracell_conc_umolL)/1000,
            intracell_conc_mmolL_sd = sd(intracell_conc_umolL)/1000)

order.of.org.df <- cul.dat2 %>% ungroup %>% 
  select(Org_Name, Org_Type, Org_Type_Specific, CultureID_short) %>%
  unique() %>% arrange(Org_Type, Org_Type_Specific)

order.of.org <- order.of.org.df$Org_Name

cul.dat2$Org_Name <- factor(cul.dat2$Org_Name, 
                                   levels = (order.of.org))


#Load up the quan enviro data-----
dat <- read_csv(quandat.file) %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)
dat2 <- dat %>%
  select(Identification, SampID,  nmolinEnviroave, nmolCave, molFractionC, RankPercent ) %>%
  left_join(meta.dat.enviro) %>%
  group_by(Station_1, Depth, Cruise, Identification) %>%
  summarise(Lat = mean(latitude),
            nmolinEnviroave_mean = mean(nmolinEnviroave, na.rm = T), 
            nmolinEnviroave_sd = sd(nmolinEnviroave, na.rm = T)) %>%
  filter(Identification %in% cmpds) 
CV.ave <- dat2 %>%
  mutate(CV = nmolinEnviroave_sd/nmolinEnviroave_mean) %>%
  group_by(Identification) %>%
  summarise(CV = mean(CV, na.rm = TRUE))
dat2 <- dat2 %>%
  left_join(CV.ave) %>%
  mutate(nmolinEnviroave_sd = ifelse(is.na(nmolinEnviroave_sd), CV*nmolinEnviroave_mean, nmolinEnviroave_sd))


#Plot up depth profiles; latitudinal gradients----
g.dpnorth <- ggplot(data = dat2 %>% filter(Cruise == "MGL1704") %>%
                      mutate(nmolinEnviroave_mean = 
                               ifelse(Identification == cmpds[1], nmolinEnviroave_mean, nmolinEnviroave_mean*transform_factor)),
                    aes(x = Depth, y = nmolinEnviroave_mean, color = Identification)) +
  geom_point()+
  geom_line() +
  scale_y_continuous(sec.axis = sec_axis(~ . / transform_factor, name = "nM Trigonelline"))+
  scale_x_reverse()+
  scale_color_manual(values = c("#1B9E77", "#B7469B"))+
  coord_flip()+
  labs(y= "nM Homarine", x = "Depth (m)") +
  theme(axis.title.x.top = element_text(size = 7, color = "#B7469B"),
        axis.title.x.bottom = element_text(size = 7, color = "#1B9E77"),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "none")

g.dpsouth <- ggplot(data = dat2 %>% filter(Cruise == "KM1513") %>%
                      mutate(nmolinEnviroave_mean = 
                               ifelse(Identification == cmpds[1], nmolinEnviroave_mean, nmolinEnviroave_mean*transform_factor),
                             nmolinEnviroave_sd = 
                               ifelse(Identification == cmpds[1], nmolinEnviroave_sd, nmolinEnviroave_sd*transform_factor)),
                    aes(x = Depth, y = nmolinEnviroave_mean, color = Identification)) +
  geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = nmolinEnviroave_mean-nmolinEnviroave_sd, 
                  ymax = nmolinEnviroave_mean+nmolinEnviroave_sd, fill = Identification),
              color = NA, alpha = 0.2)+
  scale_y_continuous(sec.axis = sec_axis(~ . / transform_factor, name = "nM Trigonelline"))+
  scale_x_reverse()+
  scale_color_manual(values = c("#1B9E77", "#B7469B"))+
  scale_fill_manual(values = c("#1B9E77", "#B7469B"))+
  coord_flip() +
  labs(y= "nM Homarine", x = "Depth (m)") +
  theme(axis.title.x.top = element_text(size = 7, color = "#B7469B"),
        axis.title.x.bottom = element_text(size = 7, color = "#1B9E77"),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "none")

g.dptransect <- ggplot(data = dat2 %>% filter(Cruise == "KOK1606") %>%
                         mutate(nmolinEnviroave_mean = 
                                  ifelse(Identification == cmpds[1], nmolinEnviroave_mean, nmolinEnviroave_mean*transform_factor),
                                nmolinEnviroave_sd = 
                                  ifelse(Identification == cmpds[1], nmolinEnviroave_sd, nmolinEnviroave_sd*transform_factor)),
                    aes(x = Lat, y = nmolinEnviroave_mean, color = Identification)) +
  geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = nmolinEnviroave_mean-nmolinEnviroave_sd, 
                  ymax = nmolinEnviroave_mean+nmolinEnviroave_sd, fill = Identification),
              color = NA, alpha = 0.2)+
  scale_y_continuous(sec.axis = sec_axis(~ ./ transform_factor, name = "nM Trigonelline")) +
  scale_color_manual(values = c("#1B9E77", "#B7469B"))+
  scale_fill_manual(values = c("#1B9E77", "#B7469B"))+
  labs(y= "nM Homarine", x = "Latitude") +
  theme(axis.title.y.right = element_text(size = 7, color = "#B7469B"),
        axis.title.y.left = element_text(size = 7, color = "#1B9E77"),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = c(0.3, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, 'lines'))

#Make chromats-----
chromat.dat2 <- chromat.dat %>%
  filter(RT < 11) %>%
  mutate(intensity = ifelse(Column == "RP", intensity*8+1E7, intensity)) %>%
  filter(RT < 5 | Column == "HILIC") %>%
  filter(RT > 4 | Column == "RP")

g.chromat <- ggplot(data =chromat.dat2, 
                    aes(y = intensity, x = RT, color = Column)) +
  geom_line() +
  scale_color_manual(values = c("grey30", "grey60"))+
  labs(y= "Intensity", x = "Retention Time (minutes)") +
  annotate("text", x = c(6.3,8.5), y = c(2.2E8,6E7), label = c("homarine","trigonelline"), 
           fontface = "italic", size = 2.5)+
  theme(axis.title = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 7))
g.chromat


#Get plot of homarine and trig in organisms------
#Munge data to get rid of orgs that we never saw these compounds in, get just species name
good.orgs <- cul.dat2 %>%
  group_by(Org_Name) %>%
  summarise(intracell_conc_mmolL_sum = sum(intracell_conc_mmolL, na.rm = FALSE)) %>%
  filter(intracell_conc_mmolL_sum>0)
cul.dat3 <- cul.dat2 %>%
  filter(Org_Name %in% good.orgs$Org_Name) %>%
  left_join(meta.dat.culture %>% select(CultureID_short, Species), by = "CultureID_short") %>%
  arrange(Org_Type_Specific)
cul.dat3$Species <- factor(cul.dat3$Species, 
                            levels = (unique(cul.dat3$Species)))

my_x_titles <- c(expression(paste(italic("Synechococcus ")~'WH7803')),
                 expression(paste(italic("Synechococcus ")~'WH8102')),
                 expression(paste(italic("Cyclotella meneghiniana"))),
                  expression(paste(italic("Pseudo-nitzschia pungens"))),
                  expression(paste(italic("Phaeodactylum tricornutum"))),
                  expression(paste(italic("Thalassiosira oceanica"))),
                 expression(paste(italic("Thalassiosira pseudonana"))),
                  expression(paste(italic("Amphidinium carterae"))),
                  expression(paste(italic("Alexandrium tamarense"))),
                 expression(paste(italic("Lingulodinium polyedra"))),
                 expression(paste(italic("Heterocapsa triquetra"))),
                  expression(paste(italic("Emiliania huxleyi ")~'CCMP2090')),
                  expression(paste(italic("Emiliania huxleyi ")~'CCMP371')),
                  expression(paste(italic("Micromonas pusilla"))))



g.cul <- ggplot(data = cul.dat3 %>%
                  mutate(intracell_conc_mmolL = 
                           ifelse(intracell_conc_mmolL == 0, NA, intracell_conc_mmolL)) %>%
                  mutate(intracell_conc_mmolL = 
                           ifelse(Identification == cmpds[1], intracell_conc_mmolL, intracell_conc_mmolL*100)),
                aes(y = Org_Name, x = intracell_conc_mmolL, fill = Identification)) +
  geom_bar(stat= "identity", position = "dodge")+
  scale_y_discrete(labels = my_x_titles)+
  scale_x_continuous(sec.axis = sec_axis(~ . / 100, name = "mM intracellular Trigonelline"), 
                     expand = c(0, 0), limits = c(0, 5E2),
                     breaks = c(0, 1E2, 2E2, 3E2,  4E2))+
  scale_fill_manual(values = c("#1B9E77", "#B7469B"))+
  labs(x= "mM intracellular Homarine") +
  theme(axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 6),
        legend.position = "none")
  
g.cul

cul.dat4 <- cul.dat3 %>%
  filter(intracell_conc_mmolL > 0) %>%
  filter(!is.na(intracell_conc_mmolL)) %>%
  unique()
   mutate(intracell_conc_mmolL = 
           ifelse(Identification == cmpds[1], intracell_conc_mmolL, intracell_conc_mmolL*100))


g.cul.2 <- ggplot(data = cul.dat4,
                  aes(y = Org_Name, x = Identification, fill = Identification, 
                      label = signif(intracell_conc_mmolL, digits = 2))) +
  geom_tile(aes(alpha = log10(intracell_conc_mmolL)))+
  geom_text(size = 2.5) +
  scale_fill_manual(values = c("#1B9E77", "#B7469B"))+
  scale_y_discrete(labels = my_x_titles,
                   expand = c(0, 0))+
  scale_x_discrete(position = "top", name="Intraceulluar concentrations (mM)",
                   expand = c(0, 0))+
  theme(axis.title.x = element_text(size = 8),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        legend.position = "none")



#Combine all the plots -----
g.combo <- plot_grid(g.dptransect, g.dpnorth, g.dpsouth, ncol = 3, rel_widths = c(2,1,1), labels = c("A", "B", "C"))

logo_file <- system.file("extdata", "logo.png", package = "cowplot")
molecule_file <- "Figures/Molecules/Homarine.pdf"
molecule_file2 <- "Figures/Molecules/Trigonelline.pdf"

g.combo2 <- ggdraw(g.combo) + 
  draw_image(molecule_file, x = 0.2, y = 1, hjust = 1, vjust = 1, width = 0.13, height = 0.2)+
  draw_image(molecule_file2, x = 0.2, y = 0.8, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
g.combo2

g.combo3 <- plot_grid(g.chromat, g.cul.2, ncol = 2, rel_widths = c(1,1.6), labels = c("D", "E"))
g.combo4 <- plot_grid(g.combo2, g.combo3, ncol = 1, labels = c("", ""))
g.combo4

save_plot("Figures/Manuscript_figures/Homarine.pdf", g.combo4, base_height = 5, base_width = 6.5, units = "in")

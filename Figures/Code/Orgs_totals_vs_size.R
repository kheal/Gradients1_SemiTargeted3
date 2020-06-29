library(tidyverse)
library(here)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot())


dat.filename <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
Meta.dat.file <- "MetaData/CultureMetaData.csv"

#Get quan data in shape to plot
dat <- read_csv(dat.filename)
meta.dat <- read_csv(Meta.dat.file)

dat.mean <- dat %>% ungroup () %>%
  group_by(Identification, ID, Org_Name, Org_Type_Specific, Org_Type) %>%
  mutate(intracell_conc_umolL = ifelse(is.na(intracell_conc_umolL), 0 , intracell_conc_umolL)) %>%
  summarise(intracell_conc_umolCL = mean(intracell_conc_umolL))

dat.total <- dat.mean %>%
  group_by(ID) %>%
  summarise(Total_mmol_L_intracell = sum(intracell_conc_umolCL)/1000) 

meta.dat.short <- meta.dat %>% 
  select(CultureID_short, Cell_size_um3, Org_Type, Org_Type_Specific) %>% 
  rename(ID = CultureID_short) %>%
  unique()

dat.total.2 <- dat.total %>%
  left_join(meta.dat.short, by = "ID")

#Get good organism order
order.of.org.df <- dat.total.2 %>% ungroup %>% 
  select(Org_Type, Org_Type_Specific, ID) %>%
  unique() %>% arrange(Org_Type, Org_Type_Specific)

order.of.org <- order.of.org.df$ID

dat.total.2$ID <- factor(dat.total.2$ID, 
                                   levels = (order.of.org))

#Plot it up as stacked bars
pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(8)[1:8])
g.total.bars <- ggplot(data = dat.total.2, aes(x = ID, 
                                               y = Total_mmol_L_intracell,
                                               fill = Org_Type_Specific)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pal)
  
g.total.bars

#Plot it up as vs cell size
pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(8)[1:8])
g.total.scatter <- ggplot(data = dat.total.2, aes(x = Cell_size_um3, 
                                               y = Total_mmol_L_intracell,
                                               fill = Org_Type_Specific, color = Org_Type_Specific)) +
  geom_point(size = 4)+
  scale_fill_manual(values = pal)+
  scale_color_manual(values = pal)+
  scale_x_log10()+
  scale_y_log10()

g.total.scatter

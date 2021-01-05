library(tidyverse)
library(here)
library(RColorBrewer)
library(nationalparkcolors)

#TO DO: make a second plot that is the total, not stacked and stretched

dat.filename <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
Meta.dat.file <- "MetaData/CultureMetaData.csv"

#Get quan data in shape to plot
dat <- read_csv(dat.filename)
meta.dat <- read_csv(Meta.dat.file)

dat.mean <- dat %>% ungroup () %>%
  group_by(Identification, ID, Org_Name, Org_Type_Specific, Org_Type) %>%
  mutate(intracell_conc_umolL = ifelse(is.na(intracell_conc_umolL), 0 , intracell_conc_umolL)) %>%
  summarise(intracell_conc_umolL = mean(intracell_conc_umolL))%>%
  filter(Identification != "DMSP")

dat.total <- dat.mean %>%
  group_by(ID) %>%
  summarise(Total_nmolCMeasured = sum(intracell_conc_umolL))

dat.mean.molefraction <- dat.mean %>% ungroup () %>%
  left_join(dat.total, by = "ID") %>%
  mutate(molfraction = intracell_conc_umolL/Total_nmolCMeasured)

#Get good organism order
order.of.org.df <- dat %>% ungroup %>% 
  select(Org_Type, Org_Type_Specific, ID) %>%
  unique() %>% arrange(Org_Type, Org_Type_Specific)

order.of.org <- order.of.org.df$ID
order.of.org.df$ID <- factor(order.of.org.df$ID, 
                             levels = (order.of.org))


#Get good compounds and set order of compounds to highlight.  This highlights the top of each, ordered by the cumulative rank-----
#Top 4 gives 19 compounds
order.of.compounds <- dat.mean %>% ungroup %>% 
  arrange(ID, desc(intracell_conc_umolL)) %>%
  group_by(ID) %>%
  mutate(ID_rank = rank(desc(intracell_conc_umolL))) %>%
  mutate(top_ten = ifelse(ID_rank < 3, 1, NA))

order.of.compounds.2 <- order.of.compounds %>%
  ungroup() %>%
  dplyr::select(ID_rank, Identification, top_ten) %>%
  group_by(Identification) %>%
  summarise(ID_rank_sum = sum(ID_rank, na.rm = TRUE),
            top_ten = sum(top_ten, na.rm = TRUE)) %>%
  filter(top_ten > 0) %>%
  arrange(desc((top_ten)))


#Get dat.mean of just the top compounds; and dat.mean of the rest----
dat.mean.highlight <-  dat.mean %>%
  filter(Identification %in% order.of.compounds.2$Identification)

dat.mean.others <-  dat.mean %>%
  filter(!Identification %in% order.of.compounds.2$Identification) %>%
  group_by(ID, Org_Name, Org_Type_Specific, Org_Type) %>%
  summarise(intracell_conc_umolL = sum(intracell_conc_umolL)) %>%
  mutate(Identification = "all others")

dat.mean.combo <- rbind(dat.mean.highlight, dat.mean.others)

dat.mean.combo$Identification = factor(dat.mean.combo$Identification, 
                                       levels = c(order.of.compounds.2$Identification, "all others")) 

dat.mean.combo$ID = factor(dat.mean.combo$ID, 
                           levels = (order.of.org))

#Plot to highlight top 20ish in the samples----
pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(16)[1:15], rep("grey", 1))

b.all <- ggplot()+
  geom_bar(stat = "identity", data = dat.mean.combo, 
           aes(x = ID, y = intracell_conc_umolL/1000, fill = Identification), color = "black", size = 0.2, alpha = 01)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = pal)+
  labs(y = "mmol metaboltes / L \n(intracellular)" )+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="top",
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6), 
        axis.text.x = element_text(size = 6),
        plot.margin = margin(.5, .5, 0, .2, "cm"))
b.all

#Color key-----
pal2 <- park_palette("Redwoods", 5)
pal3 <- c("deepskyblue4", pal2[5], pal2[2:4], pal2[1])
g.tile.org.key <- ggplot(data = order.of.org.df, 
                         aes(x = ID, y = 1, fill = Org_Type), colour = "black") +
  geom_tile() +
  geom_text(aes(label = ID), angle = -90, size = 1.7, fontface = "italic")+
  scale_fill_manual(values = pal3)+
  xlab("Organisms (colored by broad classification)")+
  theme(axis.line.y = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = -4, r = 0, b = 0, l = 0),size = 7),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none",
        plot.margin = margin(-.4, 0, .5, -.5, "cm"))

g.tile.org.key

b.all.both <- plot_grid(b.all, g.tile.org.key, ncol = 1, rel_heights = c(1, 0.1), align = "v")
b.all.both

save_plot("Figures/Presentation_figures/stackedbar_org_nmolCmpTotal.pdf", b.all.both, base_height = 4, base_width = 6.5)

#Plot the total amount
b.all.sum <- ggplot()+
  geom_bar(stat = "identity", data = dat.mean.combo, 
           aes(x = ID, y = intracell_conc_umolL/1000, fill = Org_Type, color = Org_Type), size = 0.2)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = pal3)+
  scale_color_manual(values = pal3)+
  labs(y = "mmol C in metaboltes / L \n(intracellular)" )+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = c(0.85, 0.7),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6), 
        axis.text.x = element_text(size = 6))
b.all.sum


library(tidyverse)
library(here)
library(RColorBrewer)

dat.filename <- "Intermediates/Culture_Intermediates/Quantified_LongDat_Cultures.csv"
Meta.dat.file <- "MetaData/CultureMetaData.csv"

#TO DO: Color top 5 of each 

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

#Get good compounds to highlight (top 5 of each )
compounds.to.highlight <- list()
for (i in 1:21) { 
  org = order.of.org[i]
  compounds.to.highlight[[i]] <- dat.mean %>% 
    filter(ID == org) %>%
    arrange(desc(intracell_conc_umolCL)) %>%
    head(2)
  }
compounds.to.highlight2 <- do.call(rbind, compounds.to.highlight) %>%
  ungroup() %>%
  select(Identification) %>% 
  unique()

#Get good compound order
order.of.compounds <- dat.mean %>% ungroup %>% 
  group_by(Identification) %>%
  summarise(intracell_conc_umolCL = mean(intracell_conc_umolCL, na.rm = TRUE)) %>%
  arrange(desc(intracell_conc_umolCL))


#Plot to highlight top 19 in the samples?
#pal <- c("#9B445D", rep("grey", 84))
pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(16)[1:16])
id.pal.df <- data.frame(compounds.to.highlight2$Identification, pal)
colnames(id.pal.df) <- c("Identification", "colors.avail")

dat.mean.molefraction.withcolors <- dat.mean.molefraction %>% 
  left_join(id.pal.df, by = "Identification") %>%
  mutate(color.to.plot = ifelse(Identification %in% compounds.to.highlight2$Identification, colors.avail, "white"))

dat.mean.molefraction.withcolors$Identification <- factor(dat.mean.molefraction.withcolors$Identification, 
                                               levels = (order.of.compounds$Identification))

pal.to.plot <- dat.mean.molefraction.withcolors %>%
  select(Identification, color.to.plot) %>%
  unique() %>%
  arrange(Identification)

b.all <- ggplot()+
  geom_bar(stat = "identity", data = dat.mean.molefraction.withcolors, 
           aes(x = ID, y = molfraction, fill = Identification), color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal.to.plot$color.to.plot)+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

b.all
ggsave("Figures/Preliminary/stackedbar_org.pdf", b.all, width = 8, height = 6, units = "in")



library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(plotly)
library(RColorBrewer)

#TO DO: filter out Osptists, this is normal and fine
#TO DO: filter out a few samples based on counts - plot as total to see if there are any that should be filtered out

#Set filenames (output of 18S_exploration)
otu.table.results.filename <- "Intermediates/18S/OTU_table.csv"
tax.table.results.filename <- "Intermediates/18S/tax_table.csv"
meta.dat.filename <- "Intermediates/18S/metadat_table.csv"

#Load data and join everything together
dat <- read_csv(otu.table.results.filename)
dat.long <- dat %>%
  pivot_longer(-taxon.id, names_to = "sample-id", values_to = "counts") 

meta.dat <- read_csv(meta.dat.filename)
dat.long.2 <- dat.long %>%
  left_join(meta.dat %>% mutate(`sample-id` = as.character(`sample-id`)), by = "sample-id")

taxon.dat <- read_csv(tax.table.results.filename)
dat.long.3 <- dat.long.2 %>%
  left_join(taxon.dat, by = "taxon.id")


#Get rid of samples not in the surface, dump the 28.13 lat sample, looks totally wonk
dat.long.4 <- dat.long.3 %>%
  filter(depth < 20) %>%
  filter(latitude != 28.13)

#Get names in a better place
dat.long.5 <- dat.long.4 %>%
  mutate(Org_plot = ifelse(Class == "Dinoflagellata", "Dinoflagellata", NA)) %>%
  mutate(Org_plot = ifelse(Phylum == "Alveolata" & (Class != "Dinoflagellata" | is.na(Class)), "non-Dino Alveolate", Org_plot)) %>%
  mutate(Org_plot = ifelse(Phylum == "Opisthokonta", "Opisthokonta", Org_plot))%>%
  mutate(Org_plot = ifelse(Phylum == "Rhizaria", "Rhizaria", Org_plot)) %>%
  mutate(Org_plot = ifelse(Order == "Bacillariophyta", "Diatom", Org_plot)) %>%
  mutate(Org_plot = ifelse(Phylum == "Stramenopiles" & (Order != "Bacillariophyta" |is.na(Order)), "non-Diatom Stramenopile", Org_plot)) %>%
  mutate(Org_plot = ifelse(Class == "Haptophyta", "Haptophyta", Org_plot)) %>%
  mutate(Org_plot = ifelse(Class == "Chlorophyta", "Chlorophyta", Org_plot)) %>%
  mutate(Org_plot = ifelse(Phylum == "Hacrobia", "Hacrobia", Org_plot)) %>%
  mutate(Org_plot = ifelse(Class == "Cryptophyta", "Cryptophyta", Org_plot)) %>%
  mutate(Org_plot = ifelse(is.na(Org_plot), "All Others", Org_plot))

dat.long.6 <- dat.long.5 %>%
  group_by(latitude, depth, `size-fraction`, Org_plot) %>%
  summarise(counts_sum = sum(counts))

#Better order for the Org_plot
dat.order <- dat.long.6 %>%
  group_by(Org_plot) %>%
  summarise(total_counts = sum(counts_sum)) %>%
  filter(Org_plot != "All Others") %>%
  arrange(desc(total_counts))

dat.long.6$Org_plot <- factor(dat.long.6$Org_plot, levels = c(dat.order$Org_plot, "All Others"))

#Make a bar plot, separated by size fraction
no.of.colors <- length(dat.order$Org_plot)
pal <- c(colorRampPalette(brewer.pal(7,"Dark2"))(no.of.colors)[1:no.of.colors], "grey")

b.all <- ggplot()+
  geom_bar(stat = "identity", data = dat.long.6, #position = "fill",
           aes(x = factor(latitude), y = counts_sum, fill = Org_plot), color = "black", size = 0.2)+
  # scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  facet_grid(rows = vars(`size-fraction`))+
  labs(y= "relative abundance", x = "latitude")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

b.all
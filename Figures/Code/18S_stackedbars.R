library(tidyverse)
library(here)
library(plotly)

#TO DO: assign names like in the HMMR plots - right now the names are a mess with all 9 levels
dat.filename <- "Intermediates/18S/pr3-dada2_level-9.csv"

#Get quan data in shape to plot
dat <- read_csv(dat.filename)
dat.long <- dat %>%
  select(-index) %>%
  pivot_longer(-(depth:`size-fraction`), names_to = "Orgs", values_to = "counts") %>%
  filter(depth < 20)
dat.long.2 <- dat.long %>%
  mutate(latitude = as.factor(latitude)) %>%
  group_by(latitude, Orgs, depth, `size-fraction`) %>%
  summarise(counts = mean(counts))%>%
  separate(Orgs, into = c("level.1", "level.2", "level.3", "level.4",
                          "level.5", "level.6", "level.7", "level.8", "level.9"), sep = ";", remove = FALSE) %>%
  filter(latitude != 28.13)

#Re-assign names that make more sense
dat.long.3 <- dat.long.2 %>%
  mutate(Org_plot = ifelse(level.3 == "Dinoflagellata", "Dinoflagellata", NA)) %>%
  mutate(Org_plot = ifelse(level.2 == "Alveolata" & level.3 != "Dinoflagellata", "non-Dino Alveolate", Org_plot)) %>%
  mutate(Org_plot = ifelse(level.2 == "Opisthokonta", "Opisthokonta", Org_plot))%>%
  mutate(Org_plot = ifelse(level.2 == "Rhizaria", "Rhizaria", Org_plot)) %>%
  mutate(Org_plot = ifelse(level.4 == "Bacillariophyta", "Diatom", Org_plot)) %>%
  mutate(Org_plot = ifelse(level.2 == "Stramenopiles" & level.4 != "Bacillariophyta", "non-Diatom Stramenopile", Org_plot)) %>%
  mutate(Org_plot = ifelse(level.3 == "Haptophyta", "Haptophyta", Org_plot)) %>%
  mutate(Org_plot = ifelse(level.3 == "Chlorophyta", "Chlorophyta", Org_plot)) %>%
  mutate(Org_plot = ifelse(level.3 == "Cryptophyta", "Cryptophyta", Org_plot)) %>%
  mutate(Org_plot = ifelse(is.na(Org_plot), "Others", Org_plot))

dat.long.3.inspect <- dat.long.3 %>%
  filter(level.2 != "Alveolata") %>%
  filter(level.2 != "Opisthokonta") %>%
  filter(level.2 != "Rhizaria") %>%
  filter(level.3 != "Haptophyta") %>%
  filter(level.2 != "Stramenopiles")

dat.long.4 <- dat.long.3 %>%
  group_by(latitude, depth, `size-fraction`, Org_plot) %>%
  summarise(counts_sum = sum(counts))

#Better order for the Org_plot
dat.order <- dat.long.4 %>%
  group_by(Org_plot) %>%
  summarise(total_counts = sum(counts_sum)) %>%
  arrange(desc(total_counts))

dat.long.4$Org_plot <- factor(dat.long.4$Org_plot, levels = dat.order$Org_plot)

#Make a bar plot, separated by size fraction
pal <- c(colorRampPalette(brewer.pal(7,"Dark2"))(15)[1:15])

b.all <- ggplot()+
  geom_bar(stat = "identity", data = dat.long.4, position = "fill",
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


ggsave("Figures/Preliminary/stackedbar_org.pdf", b.all, width = 8, height = 6, units = "in")


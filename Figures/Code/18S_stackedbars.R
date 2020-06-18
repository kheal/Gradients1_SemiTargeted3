library(tidyverse)
library(here)
library(plotly)

#TO DO: assign names like in the HMMR plots
dat.filename <- "Intermediates/18S/pr3-dada2_level-9.csv"

#Get quan data in shape to plot
dat <- read_csv(dat.filename)
dat.long <- dat %>%
  select(-index) %>%
  pivot_longer(-(depth:`size-fraction`), names_to = "Orgs", values_to = "rel_abun") %>%
  filter(depth < 20)
dat.long.2 <- dat.long %>%
  mutate(latitude = as.factor(latitude)) %>%
  group_by(latitude, Orgs, depth, `size-fraction`) %>%
  summarise(rel_abun = mean(rel_abun))%>%
  separate(Orgs, into = c("level 1", "level 2", "level 3", "level 4",
                          "level 5", "level 6", "level 7", "level 8", "level 9"), sep = ";", remove = FALSE) 

#Make a bar plot, separated by 
b.all <- ggplot()+
  geom_bar(stat = "identity", data = dat.long.2, position = "fill",
           aes(x = factor(latitude), y = rel_abun, fill = Orgs), color = "black", size = 0.2)+
 # scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
 # scale_fill_manual(values = pal.to.plot$color.to.plot)+
  facet_grid(rows = vars(`size-fraction`))+
  labs(y= "mol fraction C")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position="bottom",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

b.all
ggsave("Figures/Preliminary/stackedbar_org.pdf", b.all, width = 8, height = 6, units = "in")


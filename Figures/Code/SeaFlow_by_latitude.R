#Based of scripts I found here: "⁨0_Manuscript▸ ⁨_Heal_GradientsChapter⁩ ▸ ⁨_WorkingManuscript⁩ ▸ ⁨Figures/⁨Sample_Maps⁩"
library(tidyverse)
library(cowplot)
library(here)

round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

#Name your inputs
dat.file <- "MetaData/SeaFlow/SeaFlow_withKOK.csv"

#Load Seaflow dat
dat <- read_csv(dat.file)

#Make a 0.1 degree running mean
dat2 <- dat %>%
  mutate(NBorSB = ifelse(X1 > 16514, "SB", "NB"),
         BinLat = round_any(latitude, 0.01)) %>%
  group_by(pop, NBorSB, BinLat) %>%
  summarise(n_count = median(n_count, na.rm = TRUE), 
            temp1 = mean(temp1, na.rm = TRUE)) %>%
  filter(temp1>10 & temp1 < 25) %>% filter(NBorSB == "NB") %>%
  filter(pop != "beads")%>%
  filter(pop != "unknown") %>%
  mutate(pop2 = ifelse(pop == "picoeuks", "Picoeukaryotes", 
                       ifelse(pop == "prochloro", "Prochlorococcus", "Synechococcus")))

g <- ggplot(data = dat2, aes(x = BinLat, y = n_count)) +
  geom_point()+
  geom_line()+
  facet_wrap(pop2 ~ ., ncol = 1, scales = "free_y") +
  labs( x = "Latitude", y = "cells per mL")+
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        legend.position = "none")
               
g



save_plot("Figures/Manuscript_figures/Seaflow_lat.pdf", g, base_height = 6, base_width = 6.5, device=cairo_pdf)


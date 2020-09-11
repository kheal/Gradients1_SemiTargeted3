library(tidyverse)
library(chron)
library(cowplot)
theme_set(theme_cowplot())

#Get your layout info here
dat_layout <- read_delim("RawOutput/DSS3_growthcurve/RawData/20181102_DSS3_CarbonSources1.txt", "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 46) %>%
  head(16)%>%
  fill(X1) %>%
  filter(!str_detect(`1`, "SPL")) %>%
  gather("Column", "ID", -X1) %>%
  mutate(Well = paste(X1, Column, sep = "")) %>%
  select(Well, ID)

#Get your first 24 hours here
dat_first24 <- read_delim("RawOutput/DSS3_growthcurve/RawData/20181102_DSS3_CarbonSources1.txt", "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 81) %>%
  head(721) %>%
  select(-`T\xb0 Read 2:600`) %>%
  separate(Time, into = c("H", "M", "S"), sep = ":") %>%
  mutate(Time_afterInoc = as.numeric(H)+ as.numeric(M)/60 + as.numeric(S)/60) 

#Get your second 24 hours here
dat_second24 <- read_delim("RawOutput/DSS3_growthcurve/RawData/20181102_DSS3_CarbonSources1.txt", "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 81) %>%
  tail(731) %>%
  head(626) %>% #Should be 721 %>%
  select(-`T\xb0 Read 2:600`) %>%
  separate(Time, into = c("H", "M", "S"), sep = ":") %>%
  mutate(Time_afterInoc = 24 + as.numeric(H)+ as.numeric(M)/60 + as.numeric(S)/60)

#Get your third 24 hours here
dat_third24 <- read_delim("RawOutput/DSS3_growthcurve/RawData/20181102_DSS3_CarbonSources1_part2.txt", "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 78) %>%
  head(721) %>% #Should be 721 %>%
  select(-`T\xb0 Read 1:600`) %>%
  separate(Time, into = c("H", "M", "S"), sep = ":") %>%
  mutate(Time_afterInoc = 48 + as.numeric(H)+ as.numeric(M)/60 + as.numeric(S)/60)

#Get your fourth 24 hours here
dat_fourth24 <- read_delim("RawOutput/DSS3_growthcurve/RawData/20181102_DSS3_CarbonSources1_part2.txt", "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 78) %>%
  tail(2177) %>%
  head(721) %>% #Should be 721 %>%
  select(-`T\xb0 Read 1:600`) %>%
  separate(Time, into = c("H", "M", "S"), sep = ":") %>%
  mutate(Time_afterInoc = 72 + as.numeric(H)+ as.numeric(M)/60 + as.numeric(S)/60)

#Connect and turn into long data
dat_all <- rbind(dat_first24, dat_second24 ) %>%
  rbind(dat_third24) %>%
  rbind(dat_fourth24) %>%
  select(Time_afterInoc, everything()) %>%
  select(-H, -M, -S)%>%
  gather(key = "Well", value = "rawOD", -Time_afterInoc) %>%
  left_join(dat_layout)

#Get rolling media
dat_med <- dat_all %>%
  mutate(averageTime = plyr::round_any(Time_afterInoc, accuracy = 2))%>%
  group_by(Well, ID, averageTime) %>%
  summarise(medOD = median(as.numeric(rawOD)))

#Subtract no inoculumn blank
dat_blanked <- dat_med %>%
  left_join(dat_med %>%
              filter(ID == "L1_noC_neg") %>%
              group_by(averageTime) %>%
              summarise(BlankmedOD = median(as.numeric(medOD)))) %>%
  mutate(medODblanked = medOD - BlankmedOD)

#Get averaged OD for each time point
dat_averaged <- dat_blanked %>%
  group_by(ID, averageTime) %>%
  summarise(AveOD = mean(medODblanked),
            sdOD = sd(medODblanked)) %>%
  filter(!is.na(AveOD))

#Plot it up for Anitra
dat_averaged_toPlot <- dat_averaged %>%
  filter(ID %in% c("L1_noC_pos", "L1_Acetate_pos", "L1_Trigonelline_pos", "L1_Homarine_pos" ))

order <- c("L1_noC_pos", "L1_Acetate_pos", "L1_Trigonelline_pos", "L1_Homarine_pos" )

dat_averaged_toPlot$ID <- factor(dat_averaged_toPlot$ID, levels = order)

g <- ggplot(dat = dat_averaged_toPlot, aes(x = averageTime, y = AveOD, fill = ID)) +
  geom_line() +
  geom_ribbon(aes(ymin = AveOD - sdOD,
                  ymax = AveOD + sdOD), alpha = 0.5) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(labels = c("no additional \ncarbon", "Acetate", "Trigonelline", "Homarine"))+
  labs(x = "Time after inoculation (hr)", 
       y = "Optical Density (OD)", 
       fill= "Carbon source") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title =  element_text(size = 7),
        legend.text = element_text(size = 6))

g

save_plot("Figures/Manuscript_figures/DSS3_growthcurve.pdf", g, base_height = 4, base_width = 6, units = "in")


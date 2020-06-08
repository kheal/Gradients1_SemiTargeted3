library(here)
library(tidyverse)
library(xtable)
library(lubridate)

samp.info <- read_csv("MetaData/CultureMetaData.csv")
samp.info.short <- samp.info %>%
  select(Species, strain, Org_Type) %>%
  group_by(Species, strain, Org_Type) %>%
  summarise(n = n()) %>% ungroup() %>%
  filter(!is.na(strain)) %>%
  rename(Group = Org_Type, Strain = strain) %>%
  select(Group, Species, Strain, n) %>%
  arrange(Species) %>%
  arrange(Group) %>%
  mutate(Species = paste0("\\textit{", Species, "}")) %>%
  rename(`\\textit{n}` = n)

table.latex <- xtable(samp.info.short, sanitize.text.function = identity)
caption(table.latex) <- "\\label{CultureSampleDescriptions}Summary of cultured phytoplankton collected and analyzed in this study. More information (including all culturing conditions) can be found in \\cite{Durham2019}."

print(table.latex, type="latex", file="Tables/Manuscript_tables/CultureSampleDescriptions.tex",  sanitize.text.function = identity, include.rownames=FALSE)

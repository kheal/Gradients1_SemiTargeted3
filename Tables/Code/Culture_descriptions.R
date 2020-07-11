library(here)
library(tidyverse)
library(xtable)
library(lubridate)

samp.info <- read_csv("MetaData/CultureMetaData.csv")
samp.info.short <- samp.info %>%
  select(Species, strain, Org_Type, CultureID_short) %>%
  group_by(Species, strain, Org_Type, CultureID_short) %>%
  summarise(n = n()) %>% ungroup() %>%
  filter(!is.na(strain)) %>%
  rename(`Broad taxon` = Org_Type, Strain = strain, `Short ID` = CultureID_short) %>%
  select(`Broad taxon`, Species, Strain, n, `Short ID`) %>%
  arrange(Species) %>%
  arrange(`Broad taxon`) %>%
  mutate(Species = paste0("\\textit{", Species, "}")) %>%
  rename(`\\textit{n}` = n)

table.latex <- xtable(samp.info.short, sanitize.text.function = identity)
caption(table.latex) <- "\\label{CultureSampleDescriptions}Summary of cultured organisms analyzed in this study. More information (including culturing conditions for all except the Archaea) can be found in \\cite{Durham2019}. More detailed information are in Table \ref{FullCultureSampleDescriptions}. Short ID is how the organism is labeled throughout the figures."

print(table.latex, type="latex", file="Tables/Manuscript_tables/CultureSampleDescriptions.tex",  sanitize.text.function = identity, include.rownames=FALSE)

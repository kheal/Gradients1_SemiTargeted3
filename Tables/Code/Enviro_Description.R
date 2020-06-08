library(tidyverse)
library(xtable)
library(lubridate)

enviro.dat <- read_csv("MetaData/Enviro_data_fromRosie.csv")
enviro.dat.short <- enviro.dat %>%
  select(-`Number of sampled locations`) %>%
  rename(`\\makecell{Latitude \\\\ ($^\\circ$N)}` = `Latitude of stations`,
         `\\makecell{SST \\\\ ($^\\circ$C)}` = SST,
         `\\makecell{N+N \\\\ ($\\mu$M)}` = `N+N`,
         `\\makecell{PO\\textsubscript{4} \\\\ ($\\mu$M)}` = `PO4`,
         `\\makecell{Chl \\\\ (mg m\\textsuperscript{-3})}` = `Chl`)
         
table.latex <- xtable(enviro.dat.short)
caption(table.latex) <- "\\label{ZoneDescriptions} Summary of physical and chemical parameters on April 2016 cruise.  Reprinted with permission from \\cite{Gradoville2020}."

print.xtable(table.latex, type="latex", file="Tables/Manuscript_tables/RegionSummary.tex",  sanitize.text.function = identity, include.rownames=FALSE)


library(tidyverse)
library(xtable)
library(lubridate)

samp.info <- read_csv("MetaData/SampInfo_wMetaData_withUTC.csv")
samp.info.short <- samp.info %>%
  filter(!str_detect(SampID, "S2C1")) %>%
  filter(!str_detect(SampID, "S5")) %>%
  select(Cruise, Station, Depth, UTC_real, latitude, longitude, Volume, SampType) %>%
  group_by(Cruise, Station, Depth, SampType) %>%
  summarise(n = n(),
            Date = as.character(mean(as_date(UTC_real))),
            Latitude = round(mean(latitude), digits = 2),
            Longitude = round(mean(longitude), digits = 2), 
            Volume = round(mean (Volume), digits = 0)) %>% ungroup() %>%
  mutate(`\\makecell{Collection \\\\ method}` = ifelse(SampType == "U", "underway", "niskin")) %>%
  rename(`Vol (L)` = Volume, 
         `Cruise ID` = Cruise,
         `Depth (m)` = Depth) %>%
  select(-Station) %>%
  arrange(`Cruise ID`, Latitude, `Depth (m)`) %>%
  rename(`\\textit{n}` = n) %>%
  select(-`\\makecell{Collection \\\\ method}`, everything(), -SampType)

samp.info.short.withzone <- samp.info.short %>%
  mutate(`\\makecell{Environmental \\\\ regime}` = ifelse(Latitude < 31, "NPSG", "NPTZ"))


table.latex <- xtable(samp.info.short.withzone)
caption(table.latex) <- "\\label{EnviroSampleDescriptions}Summary of samples collected and analyzed in this study."

print.xtable(table.latex, type="latex", file="Tables/Manuscript_tables/EnviroSampleDescription.tex",  sanitize.text.function = identity, include.rownames=FALSE)


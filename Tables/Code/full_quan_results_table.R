# %TO DO: make supplemental table 
# \begin{table}[ht]
# \centering
# \caption{\label{FullQuanTable} Full Quantified results - each sample and compound, not a summary.} 
# \end{table}

library(tidyverse)

#Read in your dat files
dat.filename <- "Intermediates/Quantified_LongDat_Enviro.csv"

#Check out data
dat <- read_csv(dat.filename)


#Make it wide
dat.wide <- dat %>%
  select(Identification, SampID, nmolinEnviroave) %>%
  spread(key = SampID, value = nmolinEnviroave)

#Write out appropriate comment
comment <- "Quantified metabolites in all environmental samples; all values in nmol L-1"

con <- file("Tables/Manuscript_tables/SuppTables/Full_Quan_Enviro_Results.csv", open="wt")
writeLines(paste(comment), con)
write.csv( dat.wide, con)
close(con)

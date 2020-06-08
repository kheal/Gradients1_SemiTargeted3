library(tidyverse)
library(cowplot)
library(here)
library(plotly)
library(Cairo)
library(limSolve)

#File names
chemtax.input.file <- "MetaData/HPLCData/Input_Chemtax_Matrix_Mackey1998.csv"
chemtax.field.file <- "MetaData/HPLCData/KOK1606_HPLC_Chemtax_input_Mackey.csv"

#Mudge and set up the input matrix
#Input is currently based off of Fujiki 2009
chemtax.input <- read_csv(chemtax.input.file)
chemtax.input.matrix <- chemtax.input %>% select(Peri:chla) %>% as.matrix()
rownames(chemtax.input.matrix) <- chemtax.input$X1
Nx <-nrow(chemtax.input.matrix) 
EE <- rep(1, Nx)
GG <- diag(nrow = Nx)
HH <- rep(0, Nx)

chemtax.field.dat <- read_csv(chemtax.field.file)
field.samps <- chemtax.field.dat %>% select(Peri:Chla)

# Solve the matrix 
output <- lsei(t(chemtax.input.matrix), field.samps[1,] %>% as.data.frame %>% as.numeric(), EE, 1, GG, HH)$X

for (i in 2:length(field.samps$Fuco)) {
samp2.answer <- lsei(t(chemtax.input.matrix), field.samps[i,] %>% as.data.frame %>% as.numeric(), EE, 1, GG, HH)$X
output <- rbind(output, samp2.answer)
}

rownames(output) <- NULL
output <- as.tibble(output)
field.output <- cbind(chemtax.field.dat %>% select(latitude:Cast), output)

field.output.long <- field.output %>%
  gather(key = "org", value = "fraction", Chlorophytes:Prochlorphytes) %>%
  filter(!is.na(latitude)) %>%
  filter(`Depth (m)` < 20)

field.output.long.summary <- field.output.long %>%
  mutate(latitude_short = as.factor(round(latitude, digits = 1))) %>%
  group_by(latitude_short, org) %>%
  summarise(fraction = mean(fraction, na.rm = T))

#Normalize this to Chla values!
g <- ggplot(field.output.long.summary, aes(x = latitude_short, y = fraction, fill = org))+
  geom_bar(stat = "identity")
  

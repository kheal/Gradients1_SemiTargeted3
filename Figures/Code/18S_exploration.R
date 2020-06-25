#Based on Ben Lamberts github at https://github.com/armbrustlab/gradients-amplicon/blob/master/R/G1-comparison.R

library("phyloseq"); packageVersion("phyloseq")
library(qiime2R) #did not install new versions of dependencies while installing
library("ggplot2"); packageVersion("ggplot2")
library("scales")

#Set file names
dada.filename <-"Intermediates/Raw_18S/gradients-amplicon/g1/pr2-dada2-final.qza"
taxon.filename <-"Intermediates/Raw_18S/gradients-amplicon/g1/g1-taxonomy-pr2-dada2.qza"
meta.filename <-"Intermediates/Raw_18S/gradients-amplicon/g1/dada2-metadata.txt"
meta.dat <- read_delim(meta.filename, "\t", escape_double = FALSE, trim_ws = TRUE)

#Get physeq data;
physeq<-qza_to_phyloseq(
  features=dada.filename,
  taxonomy=taxon.filename,
  metadata = meta.filename
)

#Munging to get the data into a workable shape for putting into the 18S_stackedbars
otu.table.results <- physeq@otu_table@.Data %>% as.data.frame()
tax.table.results <- physeq@tax_table@.Data %>% as.data.frame()

#Writing data that can easily be plotted in 18S_stackedbars.R
write_csv(otu.table.results, "Intermediates/18S/OTU_table.csv")
write_csv(tax.table.results, "Intermediates/18S/tax_table.csv")
write_csv(meta.dat, "Intermediates/18S/metadat_table.csv")


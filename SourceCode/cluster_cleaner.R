#TO DO: clean up cluster letters - make Org cluster 
library(tidyverse)

#Read in your dat files
Org.clusterfilename <- "Intermediates/organs_wide_stand_withclusters_tofix.csv"
MGL.clusterfilename <- "Intermediates/MGL_wide_stand_withclusters_tofix.csv"

#Rename Org cluster names in MF file 
old.clusters <- c("a", "b", "c", "d", "e", "f", "g", "h")
new.clusters <- c("a", "c", "b", "g", "d", "e", "f", "h")
org.cluster.namer <-  data.frame(old.clusters,new.clusters)
names(org.cluster.namer) <- c("cluster_letters", "cluster_letters_new")
Org.dat.og <- read_csv(Org.clusterfilename)
Org.dat <- Org.dat.og%>%
  left_join(org.cluster.namer, by = "cluster_letters") %>%
  rename(cluster_letters_og = cluster_letters) %>%
  rename(cluster_letters = cluster_letters_new) %>%
  select(-cluster_letters_og)
write_csv(Org.dat, "Intermediates/organs_wide_stand_withclusters.csv")


#Rename MGL cluster names in MF file 
old.clusters <- c("a", "b", "c", "d", "e", "f", "g")
new.clusters <- c("a", "b", "c", "f", "d", "e", "f")
org.cluster.namer <-  data.frame(old.clusters,new.clusters)
names(org.cluster.namer) <- c("cluster_letters", "cluster_letters_new")
MGL.dat.og <- read_csv(MGL.clusterfilename)
MGL.dat <- MGL.dat.og%>%
  left_join(org.cluster.namer, by = "cluster_letters") %>%
  rename(cluster_letters_og = cluster_letters) %>%
  rename(cluster_letters = cluster_letters_new) %>%
  select(-cluster_letters_og)
write_csv(MGL.dat, "Intermediates/MGL_wide_stand_withclusters.csv")


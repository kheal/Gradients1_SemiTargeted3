library(here)
library(tidyverse)

#Load data----
MGL.dat.file <- "Intermediates/MGL_wide_stand_withclusters.csv"
KM.dat.file <- "Intermediates/KM_wide_stand_withclusters.csv"
KOK.dat.file <- "Intermediates/KOK_wide_stand_withclusters.csv"
Org.dat.file <- "Intermediates/organs_wide_stand_withclusters.csv"


#Get KOK cluster info----
KOK.dat <- read_csv(KOK.dat.file)

KOK.dat <- KOK.dat %>%
  mutate(KOK_clusters = paste0("KOK_", cluster_letters)) %>% 
  select(MassFeature_Column, KOK_clusters) %>%
  unique()

#Get Org cluster info----
Org.dat <- read_csv(Org.dat.file) %>%
  select(MassFeature_Column, cluster_letters) %>%
  unique() %>%
  right_join(KOK.dat, by = "MassFeature_Column") %>%
  mutate(Org_clusters = paste0("Org_", cluster_letters)) %>% 
  select(-KOK_clusters)
  

#Get MGL cluster info----
MGL.dat <- read_csv(MGL.dat.file) %>%
  select(MassFeature_Column, cluster_letters) %>%
  unique() %>%
  right_join(KOK.dat, by = "MassFeature_Column") %>%
  mutate(MGL_clusters = paste0("MGL_", cluster_letters)) %>% 
  select(-KOK_clusters)

#Get KM cluster info----
KM.dat <- read_csv(KM.dat.file)  %>%
  select(MassFeature_Column, cluster_letters) %>%
  unique() %>%
  right_join(KOK.dat, by = "MassFeature_Column") %>%
  mutate(KM_clusters = paste0("KM_", cluster_letters)) %>% 
  select(-KOK_clusters)


#Bootstrap the overlap between KOK and Org--------
#Put together KOK and Org 
MF.clusters <- KOK.dat %>% left_join(Org.dat, by = "MassFeature_Column") %>%
  mutate(MF.short = row_number()) 

MF.clusters.summary <- MF.clusters %>%
  mutate(cluster_overlap = paste0(KOK_clusters, "&", Org_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(count = n())

#Set up DF with Org clusters
Org.clusters <- MF.clusters %>%
  select(MF.short, Org_clusters) %>%
  rename(MF.orgs = MF.short)

i = 1
varname <- paste0("count.", i)

MF.bootstrap<- MF.clusters %>%
  select(MF.short, KOK_clusters) %>%
  mutate(MF.orgs = sample(MF.short, replace = FALSE)) %>%
  left_join(Org.clusters, by = "MF.orgs") %>%
  mutate(cluster_overlap = paste0(KOK_clusters, "&", Org_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(!!varname := n())

#Loop it for 1000 bootstraps
for (i in 2:999) {
varname <- paste0("count.", i)
MF.bootstrap.sub<- MF.clusters %>%
  select(MF.short, KOK_clusters) %>%
  mutate(MF.orgs = sample(MF.short, replace = FALSE)) %>%
  left_join(Org.clusters, by = "MF.orgs") %>%
  mutate(cluster_overlap = paste0(KOK_clusters, "&", Org_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(!!varname := n())
MF.bootstrap <- full_join(MF.bootstrap, MF.bootstrap.sub, by = "cluster_overlap") }

#Organize and save out results
MF.pvalues <- MF.bootstrap %>%
  gather(key = "bootstrap", value = "count", -cluster_overlap) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_overlap) %>%
  summarise(pval.01 = max(count),
            pval.05 = sort(count, decreasing = TRUE)[49],
            pval.1 = sort(count, decreasing = TRUE)[99]) %>%
  left_join(MF.clusters.summary, by = "cluster_overlap")
MF.bootrapped.pvals <- MF.pvalues %>%
  mutate(pval = ifelse(count > pval.01, 0.01, "not sig")) %>%
  mutate(pval = ifelse(count > pval.05 & pval == "not sig", 0.05, pval)) %>%
  mutate(pval = ifelse(count > pval.1 & pval == "not sig", 0.1, pval))

write_csv(MF.bootrapped.pvals, "Intermediates/BootstrapResults/KOKvsOrg.csv")




#Bootstrap the overlap between KOK and MGL--------
#Put together KOK and MGL 
MF.clusters <- KOK.dat %>% left_join(MGL.dat) %>%
  mutate(MF.short = row_number()) 

MF.clusters.summary <- MF.clusters %>%
  mutate(cluster_overlap = paste0(KOK_clusters, "&", MGL_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(count = n())

#Set up DF with MGL clusters
MGL.clusters <- MF.clusters %>%
  select(MF.short, MGL_clusters) %>%
  rename(MF.MGLs = MF.short)

i = 1
varname <- paste0("count.", i)

MF.bootstrap<- MF.clusters %>%
  select(MF.short, KOK_clusters) %>%
  mutate(MF.MGLs = sample(MF.short, replace = FALSE)) %>%
  left_join(MGL.clusters, by = "MF.MGLs") %>%
  mutate(cluster_overlap = paste0(KOK_clusters, "&", MGL_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(!!varname := n())

#Loop it for 1000 bootstraps
for (i in 2:999) {
  varname <- paste0("count.", i)
  MF.bootstrap.sub<- MF.clusters %>%
    select(MF.short, KOK_clusters) %>%
    mutate(MF.MGLs = sample(MF.short, replace = FALSE)) %>%
    left_join(MGL.clusters, by = "MF.MGLs") %>%
    mutate(cluster_overlap = paste0(KOK_clusters, "&", MGL_clusters)) %>%
    group_by(cluster_overlap) %>%
    summarise(!!varname := n())
  MF.bootstrap <- full_join(MF.bootstrap, MF.bootstrap.sub, by = "cluster_overlap") }

#MGLanize and save out results
MF.pvalues <- MF.bootstrap %>%
  gather(key = "bootstrap", value = "count", -cluster_overlap) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_overlap) %>%
  summarise(pval.01 = max(count),
            pval.05 = sort(count, decreasing = TRUE)[49],
            pval.1 = sort(count, decreasing = TRUE)[99]) %>%
  left_join(MF.clusters.summary, by = "cluster_overlap")
MF.bootrapped.pvals <- MF.pvalues %>%
  mutate(pval = ifelse(count > pval.01, 0.01, "not sig")) %>%
  mutate(pval = ifelse(count > pval.05 & pval == "not sig", 0.05, pval)) %>%
  mutate(pval = ifelse(count > pval.1 & pval == "not sig", 0.1, pval))
write_csv(MF.bootrapped.pvals, "Intermediates/BootstrapResults/KOKvsMGL.csv")





#Bootstrap the overlap between KOK and KM--------
#Put together KOK and KM 
MF.clusters <- KOK.dat %>% left_join(KM.dat) %>%
  mutate(MF.short = row_number()) 

MF.clusters.summary <- MF.clusters %>%
  mutate(cluster_overlap = paste0(KOK_clusters, "&", KM_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(count = n())

#Set up DF with KM clusters
KM.clusters <- MF.clusters %>%
  select(MF.short, KM_clusters) %>%
  rename(MF.KMs = MF.short)

i = 1
varname <- paste0("count.", i)

MF.bootstrap<- MF.clusters %>%
  select(MF.short, KOK_clusters) %>%
  mutate(MF.KMs = sample(MF.short, replace = FALSE)) %>%
  left_join(KM.clusters, by = "MF.KMs") %>%
  mutate(cluster_overlap = paste0(KOK_clusters, "&", KM_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(!!varname := n())

#Loop it for 1000 bootstraps
for (i in 2:999) {
  varname <- paste0("count.", i)
  MF.bootstrap.sub<- MF.clusters %>%
    select(MF.short, KOK_clusters) %>%
    mutate(MF.KMs = sample(MF.short, replace = FALSE)) %>%
    left_join(KM.clusters, by = "MF.KMs") %>%
    mutate(cluster_overlap = paste0(KOK_clusters, "&", KM_clusters)) %>%
    group_by(cluster_overlap) %>%
    summarise(!!varname := n())
  MF.bootstrap <- full_join(MF.bootstrap, MF.bootstrap.sub, by = "cluster_overlap") }

#KManize and save out results
MF.pvalues <- MF.bootstrap %>%
  gather(key = "bootstrap", value = "count", -cluster_overlap) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_overlap) %>%
  summarise(pval.01 = max(count),
            pval.05 = sort(count, decreasing = TRUE)[49],
            pval.1 = sort(count, decreasing = TRUE)[99]) %>%
  left_join(MF.clusters.summary, by = "cluster_overlap")
MF.bootrapped.pvals <- MF.pvalues %>%
  mutate(pval = ifelse(count > pval.01, 0.01, "not sig")) %>%
  mutate(pval = ifelse(count > pval.05 & pval == "not sig", 0.05, pval)) %>%
  mutate(pval = ifelse(count > pval.1 & pval == "not sig", 0.1, pval))
write_csv(MF.bootrapped.pvals, "Intermediates/BootstrapResults/KOKvsKM.csv")






#Bootstrap the overlap between KM and MGL--------
#Put together KM and MGL 
MF.clusters <- KM.dat %>% left_join(MGL.dat, by = "MassFeature_Column") %>%
  mutate(MF.short = row_number()) 

MF.clusters.summary <- MF.clusters %>%
  mutate(cluster_overlap = paste0(KM_clusters, "&", MGL_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(count = n())

#Set up DF with MGL clusters
MGL.clusters <- MF.clusters %>%
  select(MF.short, MGL_clusters) %>%
  rename(MF.MGLs = MF.short)

i = 1
varname <- paste0("count.", i)

MF.bootstrap<- MF.clusters %>%
  select(MF.short, KM_clusters) %>%
  mutate(MF.MGLs = sample(MF.short, replace = FALSE)) %>%
  left_join(MGL.clusters, by = "MF.MGLs") %>%
  mutate(cluster_overlap = paste0(KM_clusters, "&", MGL_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(!!varname := n())

#Loop it for 1000 bootstraps
for (i in 2:999) {
  varname <- paste0("count.", i)
  MF.bootstrap.sub<- MF.clusters %>%
    select(MF.short, KM_clusters) %>%
    mutate(MF.MGLs = sample(MF.short, replace = FALSE)) %>%
    left_join(MGL.clusters, by = "MF.MGLs") %>%
    mutate(cluster_overlap = paste0(KM_clusters, "&", MGL_clusters)) %>%
    group_by(cluster_overlap) %>%
    summarise(!!varname := n())
  MF.bootstrap <- full_join(MF.bootstrap, MF.bootstrap.sub, by = "cluster_overlap") }

#MGLanize and save out results
MF.pvalues <- MF.bootstrap %>%
  gather(key = "bootstrap", value = "count", -cluster_overlap) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_overlap) %>%
  summarise(pval.01 = max(count),
            pval.05 = sort(count, decreasing = TRUE)[49],
            pval.1 = sort(count, decreasing = TRUE)[99]) %>%
  left_join(MF.clusters.summary, by = "cluster_overlap")
MF.bootrapped.pvals <- MF.pvalues %>%
  mutate(pval = ifelse(count > pval.01, 0.01, "not sig")) %>%
  mutate(pval = ifelse(count > pval.05 & pval == "not sig", 0.05, pval)) %>%
  mutate(pval = ifelse(count > pval.1 & pval == "not sig", 0.1, pval))
write_csv(MF.bootrapped.pvals, "Intermediates/BootstrapResults/KMvsMGL.csv")





#Bootstrap the overlap between KOK and Org--------
#Put together KM and Org 
MF.clusters <- KM.dat %>% left_join(Org.dat, by = "MassFeature_Column") %>%
  mutate(MF.short = row_number()) 

MF.clusters.summary <- MF.clusters %>%
  mutate(cluster_overlap = paste0(KM_clusters, "&", Org_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(count = n())

#Set up DF with Org clusters
Org.clusters <- MF.clusters %>%
  select(MF.short, Org_clusters) %>%
  rename(MF.orgs = MF.short)

i = 1
varname <- paste0("count.", i)

MF.bootstrap<- MF.clusters %>%
  select(MF.short, KM_clusters) %>%
  mutate(MF.orgs = sample(MF.short, replace = FALSE)) %>%
  left_join(Org.clusters, by = "MF.orgs") %>%
  mutate(cluster_overlap = paste0(KM_clusters, "&", Org_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(!!varname := n())

#Loop it for 1000 bootstraps
for (i in 2:999) {
  varname <- paste0("count.", i)
  MF.bootstrap.sub<- MF.clusters %>%
    select(MF.short, KM_clusters) %>%
    mutate(MF.orgs = sample(MF.short, replace = FALSE)) %>%
    left_join(Org.clusters, by = "MF.orgs") %>%
    mutate(cluster_overlap = paste0(KM_clusters, "&", Org_clusters)) %>%
    group_by(cluster_overlap) %>%
    summarise(!!varname := n())
  MF.bootstrap <- full_join(MF.bootstrap, MF.bootstrap.sub, by = "cluster_overlap") }

#Organize and save out results
MF.pvalues <- MF.bootstrap %>%
  gather(key = "bootstrap", value = "count", -cluster_overlap) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_overlap) %>%
  summarise(pval.01 = max(count),
            pval.05 = sort(count, decreasing = TRUE)[49],
            pval.1 = sort(count, decreasing = TRUE)[99]) %>%
  left_join(MF.clusters.summary, by = "cluster_overlap")
MF.bootrapped.pvals <- MF.pvalues %>%
  mutate(pval = ifelse(count > pval.01, 0.01, "not sig")) %>%
  mutate(pval = ifelse(count > pval.05 & pval == "not sig", 0.05, pval)) %>%
  mutate(pval = ifelse(count > pval.1 & pval == "not sig", 0.1, pval))

write_csv(MF.bootrapped.pvals, "Intermediates/BootstrapResults/KMvsOrg.csv")




#Bootstrap the overlap between MGL and Org--------
#Put together MGL and Org 
MF.clusters <- MGL.dat %>% left_join(Org.dat, by = "MassFeature_Column") %>%
  mutate(MF.short = row_number()) 

MF.clusters.summary <- MF.clusters %>%
  mutate(cluster_overlap = paste0(MGL_clusters, "&", Org_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(count = n())

#Set up DF with Org clusters
Org.clusters <- MF.clusters %>%
  select(MF.short, Org_clusters) %>%
  rename(MF.orgs = MF.short)

i = 1
varname <- paste0("count.", i)

MF.bootstrap<- MF.clusters %>%
  select(MF.short, MGL_clusters) %>%
  mutate(MF.orgs = sample(MF.short, replace = FALSE)) %>%
  left_join(Org.clusters, by = "MF.orgs") %>%
  mutate(cluster_overlap = paste0(MGL_clusters, "&", Org_clusters)) %>%
  group_by(cluster_overlap) %>%
  summarise(!!varname := n())

#Loop it for 1000 bootstraps
for (i in 2:999) {
  varname <- paste0("count.", i)
  MF.bootstrap.sub<- MF.clusters %>%
    select(MF.short, MGL_clusters) %>%
    mutate(MF.orgs = sample(MF.short, replace = FALSE)) %>%
    left_join(Org.clusters, by = "MF.orgs") %>%
    mutate(cluster_overlap = paste0(MGL_clusters, "&", Org_clusters)) %>%
    group_by(cluster_overlap) %>%
    summarise(!!varname := n())
  MF.bootstrap <- full_join(MF.bootstrap, MF.bootstrap.sub, by = "cluster_overlap") }

#Organize and save out results
MF.pvalues <- MF.bootstrap %>%
  gather(key = "bootstrap", value = "count", -cluster_overlap) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_overlap) %>%
  summarise(pval.01 = max(count),
            pval.05 = sort(count, decreasing = TRUE)[49],
            pval.1 = sort(count, decreasing = TRUE)[99]) %>%
  left_join(MF.clusters.summary, by = "cluster_overlap")
MF.bootrapped.pvals <- MF.pvalues %>%
  mutate(pval = ifelse(count > pval.01, 0.01, "not sig")) %>%
  mutate(pval = ifelse(count > pval.05 & pval == "not sig", 0.05, pval)) %>%
  mutate(pval = ifelse(count > pval.1 & pval == "not sig", 0.1, pval))

write_csv(MF.bootrapped.pvals, "Intermediates/BootstrapResults/MGLvsOrg.csv")




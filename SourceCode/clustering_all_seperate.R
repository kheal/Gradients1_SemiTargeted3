#This code does the clustering and writes out "Intermediates/MGL_wide_stand_withclusters.csv" files (separate) for the three cruises and for the organisms
library(tidyverse)
source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(scales)
library(cowplot)
theme_set(theme_cowplot())
library(here)

#Name your inputs
meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
numberofclusters.KOK <- 8 
numberofclusters.KM <- 4  #8 or 9 are also good
numberofclusters.MGL <- 7   
numberofclusters.org <- 7
combo.wide.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"

#Mudge data to make it usable-----
#Get names of the environmental samples
combo.wide <- read_csv(combo.wide.file)
meta.dat <- read_csv((meta.dat.file)) %>% 
  filter(SampID %in% colnames(combo.wide))
enviro.names <- meta.dat$SampID

#Get number of MFs in each HILICPos, HILICNeg, RP
combo.wide.info <- combo.wide %>% select(MassFeature_Column, Column, z) %>%
  group_by(Column, z) %>% summarise(count=n())


# Make the combo.file into a matrix
combo.wide.matrix<- combo.wide %>% select(c(enviro.names)) %>% as.data.frame()
row.names(combo.wide.matrix) <- combo.wide$MassFeature_Column

#Standardize by z scores ### for each of the cruises separately----
#Make all NAs into 0s and get rid of MFs that are all 0s, then standardize
##KM cruise first
KM.wide.matrix <- combo.wide.matrix %>%
  select(contains("KM1513"))
KM.wide.matrix[is.na(KM.wide.matrix)] <- 0
KM.wide.matrix$rowSum <- rowSums(KM.wide.matrix)
KM.wide.matrix.2 <- KM.wide.matrix %>% mutate(names = row.names(KM.wide.matrix)) %>%
  filter(!rowSum == 0) %>%
  select(-rowSum)
KM.names <-KM.wide.matrix.2$names
KM.wide.matrix.2 <- KM.wide.matrix.2 %>% select(-names)
datwidestd.noNA.KM <- data.stand((KM.wide.matrix.2), method='total', margin='row', plot=F)
row.names(datwidestd.noNA.KM) <- KM.names


##MGL cruise second
MGL.wide.matrix <- combo.wide.matrix %>%
  select(contains("S7"))
MGL.wide.matrix[is.na(MGL.wide.matrix)] <- 0
MGL.wide.matrix$rowSum <- rowSums(MGL.wide.matrix)
MGL.wide.matrix.2 <- MGL.wide.matrix %>% mutate(names = row.names(MGL.wide.matrix)) %>%
  filter(!rowSum == 0) %>%
  select(-rowSum)
MGL.names <- MGL.wide.matrix.2$names
MGL.wide.matrix.2 <- MGL.wide.matrix.2 %>% select(-names)
datwidestd.noNA.MGL <- data.stand((MGL.wide.matrix.2), method='total', margin='row', plot=F)
rownames(datwidestd.noNA.MGL) <- MGL.names


#KOK cruise third
KOK.wide.matrix <- combo.wide.matrix %>%
  select(-contains("S7")) %>%
  select(-contains("KM1513"))
KOK.wide.matrix[is.na(KOK.wide.matrix)] <- 0
KOK.wide.matrix$rowSum <- rowSums(KOK.wide.matrix)
KOK.wide.matrix <- KOK.wide.matrix %>% filter(!rowSum == 0) %>% select(-rowSum)
datwidestd.noNA.KOK <- data.stand((KOK.wide.matrix), method='total', margin='row', plot=F)
rownames(datwidestd.noNA.KOK) <- rownames(combo.wide.matrix)

#Cluster for just the KOK cruise---------
# Try to make an ordered similarity matrix in order to pick out groups from it.
datwidestd.euclid.KOK <- vegdist(datwidestd.noNA.KOK, method="euclidean") #This one is faster

# This step  gives some hints to how many clusters we should use - higher average silhoutte width the better!! Looks like for KOK either 5 or 10 clusters is best.  I'll go with 5 for now.
scree.dat.KOK <- nhclus.scree(datwidestd.euclid.KOK, max.k=50) %>% as.data.frame()
scree.plot.KOK <- ggplot(dat =scree.dat.KOK, aes(x = `no. clusters`, y = `ave silhouette width`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.KOK) 
scree.plot.2.KOK <- ggplot(dat =scree.dat.KOK, aes(x = `no. clusters`, y = `sum within-clus diss`)) +
    geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.KOK) 
scree.plot.2.KOK.combined <- plot_grid(scree.plot.KOK, scree.plot.2.KOK, ncol=1)
save_plot("Intermediates/KMeans_reports/KOK_Kmeans.pdf", scree.plot.2.KOK.combined)

#Cluster with perscribed number of clusters (set as an input)
datclu.clara.KOK <- clara(datwidestd.noNA.KOK, k = numberofclusters.KOK, metric = "euclidean", samples = 1000, sampsize = nrow(datwidestd.noNA.KOK)) 

#Add cluster info to the wide dataframe info
clusterIDs.KOK <- datclu.clara.KOK$clustering %>% as.data.frame() %>%
  rename(cluster = ".") 
clusterIDs.KOK <- clusterIDs.KOK %>%
  mutate(MassFeature_Column = rownames(clusterIDs.KOK)) 
datwidestd.noNA.KOK.wclusters <- datwidestd.noNA.KOK %>%
  mutate(MassFeature_Column = rownames(datwidestd.noNA.KOK))%>% 
  gather(., key = SampID, value = std_area, -MassFeature_Column) %>%
  mutate(SampID = SampID %>% str_replace("\\.", "\\-")) %>%
  left_join(clusterIDs.KOK) %>%
  left_join(meta.dat %>% select(SampID, latitude))

cluster <- c(1, 2, 4 ,3, 5, 6, 7, 8)
cluster_letters <- c("a", "b", "c", "d", "e", "f", "g", "h") 
KOK.mode.matcher <- cbind(cluster, cluster_letters) %>% as_tibble() %>% mutate(cluster = as.numeric(cluster))
datwidestd.noNA.KOK.wclusters <- left_join(datwidestd.noNA.KOK.wclusters, KOK.mode.matcher, by = "cluster") 

write_csv(datwidestd.noNA.KOK.wclusters, "Intermediates/KOK_wide_stand_withclusters.csv")

#Cluster for just the KM cruise---------
#Euclidean distance matrix
datwidestd.euclid.KM <- vegdist(datwidestd.noNA.KM, method="euclidean") #This one is faster

# This step  gives some hints to how many clusters we should use - higher average silhoutte width the better!! 
scree.dat.KM <- nhclus.scree(datwidestd.euclid.KM, max.k=50) %>% as.data.frame()
scree.plot.KM <- ggplot(dat =scree.dat.KM, aes(x = `no. clusters`, y = `ave silhouette width`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.KM) 
scree.plot.2.KM <- ggplot(dat =scree.dat.KM, aes(x = `no. clusters`, y = `sum within-clus diss`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.KM) 
scree.plot.2.KM.combined <- plot_grid(scree.plot.KM, scree.plot.2.KM, ncol=1)

save_plot("Intermediates/KMeans_reports/KM_Kmeans.pdf", scree.plot.2.KM.combined)

#Cluster with perscribed number of clusters (set as an input)
datclu.clara.KM <- clara(datwidestd.noNA.KM, k = numberofclusters.KM, metric = "euclidean", samples = 100, sampsize = nrow(datwidestd.noNA.KM)) 

#Make a heatmap
clusterIDs.KM <- datclu.clara.KM$clustering %>% as.data.frame() %>%
  rename(cluster = ".") 
clusterIDs.KM <- clusterIDs.KM %>%
  mutate(MassFeature_Column = rownames(clusterIDs.KM)) 

#Add cluster info to the wide dataframe info
datwidestd.noNA.KM.wclusters <- datwidestd.noNA.KM %>%
  mutate(MassFeature_Column = rownames(datwidestd.noNA.KM))%>% 
  gather(., key = SampID, value = std_area, -MassFeature_Column) %>%
  mutate(SampID = SampID %>% str_replace("\\.", "\\-")) %>%
  left_join(clusterIDs.KM) %>%
  left_join(meta.dat %>% select(SampID, latitude, Depth))

cluster <- c(1,2,3,4)
cluster_letters <- c("a", "b", "c", "d") 
KM.mode.matcher <- cbind(cluster, cluster_letters) %>% as.tibble() %>% mutate(cluster = as.numeric(cluster))
datwidestd.noNA.KM.wclusters <- left_join(datwidestd.noNA.KM.wclusters, KM.mode.matcher)

write_csv(datwidestd.noNA.KM.wclusters, "Intermediates/KM_wide_stand_withclusters.csv")

#Cluster for just the MGL cruise---------
# Try to make an ordered similarity matrix in order to pick out groups from it.
datwidestd.euclid.MGL <- vegdist(datwidestd.noNA.MGL, method="euclidean") #This one is faster

# This step  gives some hints to how many clusters we should use - higher average silhoutte width the better!! Looks like for KOK either 5 or 10 clusters is best.  I'll go with 5 for now.
scree.dat.MGL <- nhclus.scree(datwidestd.euclid.MGL, max.k=50) %>% as.data.frame()
scree.plot.MGL <- ggplot(dat =scree.dat.MGL, aes(x = `no. clusters`, y = `ave silhouette width`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.MGL) 
scree.plot.2.MGL <- ggplot(dat =scree.dat.MGL, aes(x = `no. clusters`, y = `sum within-clus diss`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.MGL) 
scree.plot.2.MGL.combined <- plot_grid(scree.plot.MGL, scree.plot.2.MGL, ncol=1)

save_plot("Intermediates/KMeans_reports/MGL_Kmeans.pdf", scree.plot.2.MGL.combined)


#Cluster with perscribed number of clusters (set as an input)
datclu.clara.MGL <- clara(datwidestd.noNA.MGL, k = numberofclusters.MGL, metric = "euclidean", samples = 100, sampsize = nrow(datwidestd.noNA.MGL)) 

#Make a heatmap
clusterIDs.MGL <- datclu.clara.MGL$clustering %>% as.data.frame() %>%
  rename(cluster = ".") 
clusterIDs.MGL <- clusterIDs.MGL %>%
  mutate(MassFeature_Column = rownames(clusterIDs.MGL)) 

#Add cluster info to the wide dataframe info
datwidestd.noNA.MGL.wclusters <- datwidestd.noNA.MGL %>%
  mutate(MassFeature_Column = rownames(datwidestd.noNA.MGL))%>% 
  gather(., key = SampID, value = std_area, -MassFeature_Column) %>%
  mutate(SampID = SampID %>% str_replace("\\.", "\\-")) %>%
  left_join(clusterIDs.MGL) %>%
  left_join(meta.dat %>% select(SampID, latitude, Depth))

cluster <- c(4,1,2,3,7,6,5)
cluster_letters <- c("a", "b", "c", "d", "e", "f", "g") 
MGL.mode.matcher <- cbind(cluster, cluster_letters) %>% as.tibble() %>% mutate(cluster = as.numeric(cluster))
datwidestd.noNA.MGL.wclusters <- left_join(datwidestd.noNA.MGL.wclusters, MGL.mode.matcher)

write_csv(datwidestd.noNA.MGL.wclusters, "Intermediates/MGL_wide_stand_withclusters.csv")

#Cluster by organism------
Org.wide.matrix<- combo.wide %>% select(-c(enviro.names)) %>% as.data.frame() %>%
  select(-Identification:-z)
row.names(Org.wide.matrix) <- Org.wide.matrix$MassFeature_Column
Org.wide.matrix <- Org.wide.matrix %>% select(-MassFeature_Column) 
Org.wide.matrix[is.na(Org.wide.matrix)] <- 0
Org.wide.matrix$rowSum <- rowSums(Org.wide.matrix)
Org.wide.matrix.2 <- Org.wide.matrix %>% mutate(names = row.names(Org.wide.matrix)) %>%
  filter(!rowSum == 0) %>%
  select(-rowSum)
Org.names <- Org.wide.matrix.2$names
Org.wide.matrix.2 <- Org.wide.matrix.2 %>% select(-names)
datwidestd.noNA.Org <- data.stand((Org.wide.matrix.2), method='max', margin='row', plot=F)
rownames(datwidestd.noNA.Org) <- Org.names

# Try to make an ordered similarity matrix in order to pick out groups from it.
datwidestd.euclid.Org <- vegdist(datwidestd.noNA.Org, method="euclidean") #This one is faster

# This step  gives some hints to how many clusters we should use - higher average silhoutte width the better!! Looks like for KOK either 5 or 10 clusters is best.  I'll go with 5 for now.
scree.dat.Org <- nhclus.scree(datwidestd.euclid.Org, max.k=50) %>% as.data.frame()
scree.plot.Org <- ggplot(dat =scree.dat.Org, aes(x = `no. clusters`, y = `ave silhouette width`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.org) 
scree.plot.2.Org <- ggplot(dat =scree.dat.Org, aes(x = `no. clusters`, y = `sum within-clus diss`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters.org) 
scree.plot.2.Org.combined <- plot_grid(scree.plot.Org, scree.plot.2.Org, ncol=1)
save_plot("Intermediates/KMeans_reports/Org_Kmeans.pdf", scree.plot.2.Org.combined)


#Cluster with perscribed number of clusters (set as an input)
datclu.clara.org <- clara(datwidestd.noNA.Org, k = numberofclusters.org, metric = "euclidean", samples = 100, sampsize = nrow(datwidestd.noNA.Org)) 

clusterIDs.org <- datclu.clara.org$clustering %>% as.data.frame() %>%
  rename(cluster = ".") 
clusterIDs.org <- clusterIDs.org %>%
  mutate(MassFeature_Column = rownames(clusterIDs.org)) 

#Add cluster info to the wide dataframe info
datwidestd.orgs.noNA.wclusters <- datwidestd.noNA.Org %>%
  mutate(MassFeature_Column = rownames(datwidestd.noNA.Org))%>% 
  gather(., key = SampID, value = std_area, -MassFeature_Column) %>%
  mutate(SampID = SampID %>% str_replace("\\.", "\\-")) %>%
  left_join(clusterIDs.org) 

cluster <- c(1,2,3,6,7,5,4)
cluster_letters <- c("a", "b", "c", "d", "e", "f", "g") 
Org.mode.matcher <- cbind(cluster, cluster_letters) %>% as_tibble() %>% mutate(cluster = as.numeric(cluster))
datwidestd.orgs.noNA.wclusters <- left_join(datwidestd.orgs.noNA.wclusters, Org.mode.matcher)


write_csv(datwidestd.orgs.noNA.wclusters, "Intermediates/organs_wide_stand_withclusters.csv")


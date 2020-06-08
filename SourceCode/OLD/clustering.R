#Name your inputs
 meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
 numberofclusters <- 10
 combo.wide.file <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
 cruise.id <- "KOK1606"

###FROM HERE DOWN WILL BE IN AN R SCRIPT AT SOME POINT
source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(scales)
library(cowplot)

#Outputs
cluster.report <- list(combo.wide.file)

combo.wide <- read_csv(combo.wide.file)

#Get names of the samples from G1
meta.dat <- read_csv(meta.dat.file) %>% 
  filter(SampID %in% colnames(combo.wide))
enviro.names <- meta.dat$SampID

# Make the combo.file into a matrix
combo.wide.matrix<- combo.wide %>% select(c(enviro.names)) %>% as.data.frame()
row.names(combo.wide.matrix) <- combo.wide$MassFeature_Column

# Standardize by z scores ### we should do this 3x for each cruise...then add them back together
KM.wide.matrix <- combo.wide.matrix %>%
  select(contains("KM1513"))
datwidestd.km <- data.stand((KM.wide.matrix), method='standardize', margin='row', plot=F)

MGL.wide.matrix <- combo.wide.matrix %>%
  select(contains("S7"))
datwidestd.MGL <- data.stand((MGL.wide.matrix), method='standardize', margin='row', plot=F)

KOK.wide.matrix <- combo.wide.matrix %>%
  select(-contains("S7")) %>%
  select(-contains("KM1513"))
datwidestd.KOK <- data.stand((KOK.wide.matrix), method='standardize', margin='row', plot=F)

datwidestd <- c(datwidestd.km, datwidestd.MGL, datwidestd.KOK ) %>% as.data.frame()
row.names(datwidestd) <- rownames(datwidestd.km)

#Need to get rid of any NANs...
datwidestd.noNA <- datwidestd %>%
  mutate(names = row.names(datwidestd)) %>%
  mutate(average = rowMeans(datwidestd)) %>%
  filter(!is.na(average)) %>% select(-average)
row.names(datwidestd.noNA) <- datwidestd.noNA$names
datwidestd.noNA <- datwidestd.noNA %>% select(-names)

# Try to make an ordered similarity matrix in order to pick out groups from it.
datwidestd.euclid <- vegdist(datwidestd.noNA, method="euclidean") #This one is faster

# This step  gives some hints to how many clusters we should use - higher average silhoutte width the better!!
scree.dat <- nhclus.scree(datwidestd.euclid, max.k=50) %>% as.data.frame()
scree.plot <- ggplot(dat =scree.dat, aes(x = `no. clusters`, y = `ave silhouette width`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters) 
scree.plot.2 <- ggplot(dat =scree.dat, aes(x = `no. clusters`, y = `sum within-clus diss`)) +
  geom_point() +  geom_line() + geom_vline(xintercept = numberofclusters) 
cluster.report[[1]] <- plot_grid(scree.plot, scree.plot.2, ncol=1)

#Cluster with perscribed number of clusters (set as an input)
datclu.clara <- clara(datwidestd.noNA, k = numberofclusters, metric = "euclidean", samples = 100, sampsize = nrow(datwidestd.noNA)) 

#Make a heatmap
clusterIDs <- datclu.clara$clustering %>% as.data.frame() %>%
  rename(cluster = ".") 
clusterIDs <- clusterIDs %>%
  mutate(MassFeature_Column = rownames(clusterIDs)) 

#Add cluster info to the wide dataframe info
datwidestd.noNA <- datwidestd.noNA %>%
  mutate(MassFeature_Column = rownames(datwidestd.noNA))%>% 
  gather(., key = SampID, value = std_area, -MassFeature_Column) %>%
  mutate(SampID = SampID %>% str_replace("\\.", "\\-")) %>%
  left_join(clusterIDs) %>%
  left_join(meta.dat %>% select(SampID, latitude))

cluster.report[[2]] <- datwidestd.noNA   
write_csv(cluster.report[[2]], "Intermediates/Enviro_longstd_withCluster.csv")


#Name your inputs
 dat.long.std.file <- "Intermediates/Enviro_longstd_withCluster.csv"
 culuture.filename <- "Intermediates/combined_wide_LogBioArea.csv"
 combo.wide.filename <- "Intermediates/WideArea_withIDinfo_withCultureLogBioArea.csv"
 meta.dat.culture.filename <- "MetaData/CultureMetaData.csv"
 numberofclusters <- 10
 meta.dat.file <- "MetaData/SampInfo_wMetaData_withUTC.csv"
 
#######

cluster.plots <- list()

 
   
#Read in your files
dat.long.std <- read_csv(dat.long.std.file)
culture.file <- read_csv(culuture.filename)
combo.wide <- read_csv(combo.wide.filename)
meta.culture.dat <- read_csv(meta.dat.culture.filename)


#Make the samples go from South to North and shallow to deep
meta.data.info <- read_csv(meta.dat.file) %>%
  select(SampID, Depth, latitude, Cruise)

orderofSamps <- dat.long.std %>% select(SampID) %>% 
  unique() %>% left_join(meta.data.info) %>%
  arrange(Depth) %>%
  arrange(latitude) %>%
  arrange(Cruise)
orderofSamps <- orderofSamps$SampID


#First lets just make a tile of culture stuff
culture.names <- colnames(culture.file)[colnames(culture.file) != "MassFeature_Column"]
culture.dat <- combo.wide %>% select(MassFeature_Column,  c(culture.names)) %>%
  gather(-MassFeature_Column, key = "CultureID", value = "LogBioVol") %>%
  left_join(combo.wide %>% select(MassFeature_Column, Identification, Confidence, mz, rt, Column, z)) %>%
  mutate(LogBioVol = ifelse(LogBioVol == 0, NA, round(as.numeric(LogBioVol), digits = 3))) %>%
  filter(!str_detect(CultureID, "P55" )) %>%
  filter(!str_detect(CultureID, "P3" ))


#Make another tile for the culture data
culture.dat.2 <- culture.dat %>% 
  left_join(dat.long.std %>% select(MassFeature_Column, cluster)) %>% 
  left_join(meta.culture.dat)  
culture.dat.3 <- culture.dat.2 %>% 
  arrange(Org_Type, Org_Type_Specific) %>% 
  select(Org_Type_Specific, CultureID) %>% unique()
culture.dat.2$CultureID <- factor(culture.dat.2$CultureID, levels = culture.dat.3$CultureID)

#Get a new order for the Mass_features based on 1) cluster and 2) number of cultures it was detected in
culture.dat.rep.num <- culture.dat.2 %>%
  filter(!is.na(LogBioVol)) %>% 
  ungroup() %>% group_by(MassFeature_Column) %>% 
  summarize(culturereps = n())

orderofIds.df <- dat.long.std %>% select(MassFeature_Column, cluster) %>% 
  unique() %>% 
  left_join(culture.dat.rep.num) %>%
  mutate(culturereps = ifelse(is.na(culturereps), 0, culturereps)) %>%
  arrange(cluster, culturereps)
orderofIds <- orderofIds.df$MassFeature_Column

dat.long.std$SampID <- factor(dat.long.std$SampID, levels = orderofSamps)
dat.long.std$MassFeature_Column <- factor(dat.long.std$MassFeature_Column, levels = rev(orderofIds))
culture.dat.2$MassFeature_Column <- factor(culture.dat.2$MassFeature_Column, levels = rev(orderofIds))

g <- ggplot(dat = dat.long.std %>%
              mutate(std_area = ifelse(std_area > 3, 3, std_area)) %>%
              mutate(std_area = ifelse(std_area < -3, -3, std_area)), aes(x = SampID, y = MassFeature_Column, fill = std_area)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  scale_fill_gradient2(low="slateblue4", mid="white", high="firebrick2",
                       midpoint=0, limits = c(-3, 3), na.value = "grey")  +
  theme(axis.text.x = element_text(angle=-60, hjust=0),
        axis.ticks = element_blank(), axis.text.y = element_blank()) +
  theme(legend.position="bottom") +
  labs(x ="Environmental samples", y = "Mass Feature")

g.cul <- ggplot(dat = culture.dat.2, aes(x = CultureID, y = MassFeature_Column, fill = LogBioVol)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  scale_fill_gradient2(low="white", high="firebrick2", na.value = "gray95")+
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 8),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())  +
  theme(legend.position="bottom")+
  labs(x ="Cultures", y = "Mass Feature")

#Combine these guys
cluster.plots[[1]] <- plot_grid(g, g.cul, ncol = 2, rel_widths = c(1, 0.5), align = "h")


#Plot again with just the medians
dat.long.std.ave <- dat.long.std %>% 
  mutate(latitude = as.factor(round(latitude, digits =1))) %>%
  group_by(latitude, MassFeature_Column, cluster) %>%
  summarize(std_area = median(std_area))

culture.dat.2.ave <- culture.dat.2 %>% 
  group_by(CultureID_short, MassFeature_Column, Org_Type, Org_Type_Specific, cluster) %>%
  mutate(LogBioVol = ifelse(is.na(LogBioVol), 0, LogBioVol)) %>%
  summarize(LogBioVol = median(LogBioVol)) %>%
  mutate(LogBioVol = ifelse(LogBioVol == 0, NA, LogBioVol))

culture.dat.rep.num.ave <- culture.dat.2.ave %>%
  filter(!is.na(LogBioVol)) %>% 
  ungroup() %>% group_by(MassFeature_Column, cluster) %>% 
  summarize(culturereps = n()/21)
culture.dat.rep.num.ave$MassFeature_Column <- factor(culture.dat.rep.num.ave$MassFeature_Column, levels = rev(orderofIds))


culture.dat.3.ave <- culture.dat.2.ave %>% ungroup() %>% arrange(Org_Type, Org_Type_Specific) %>% 
  select(CultureID_short) %>% unique()
culture.dat.2.ave$CultureID_short <- factor(culture.dat.2.ave$CultureID_short, levels = culture.dat.3.ave$CultureID_short)



g2 <- ggplot(dat = dat.long.std.ave %>% 
               mutate(std_area = ifelse(std_area > 2.5, 2.5, std_area)) %>%
               mutate(std_area = ifelse(std_area < -2.5, -2.5, std_area)), aes(x = latitude, y = MassFeature_Column, fill = std_area)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  scale_fill_gradient2(low="slateblue4", mid="white", high="firebrick2",
                       midpoint=0, limits = c(-3, 3), na.value = "grey")  +
  theme(axis.text.x = element_text(angle=-60, hjust=0),
        axis.ticks = element_blank(), axis.text.y = element_blank(), 
        strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(legend.position="bottom") +
  labs(x ="Latitude", y = "Mass Feature", fill = "Standardized \n peak area")

g.cul2 <- ggplot(dat = culture.dat.2.ave, aes(x = CultureID_short, y = MassFeature_Column, fill = LogBioVol)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  scale_fill_gradient2(low="white", high="seagreen4", na.value = "gray95")+
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 8),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())  +
  theme(legend.position="bottom") +
  labs(x ="Culture ID", y = "Mass Feature", fill = "Log(Peak area/biovolume)")


g.culseen <- ggplot(dat = culture.dat.rep.num.ave, aes(x = 1, y = MassFeature_Column, fill = culturereps)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales="free_y", space="free_y")+
  scale_fill_gradient2(low="white", high="black")+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(), 
        strip.text.y = element_blank())  +
  theme(legend.position="bottom") 


#Combine these guys (meidans of station)
cluster.plots[[2]]  <- plot_grid(g2, g.culseen, g.cul2, ncol = 3, rel_widths = c(1, 0.05, 0.5), align = "h")



#Make plot that is average mode - what the average compound looks like in each of these clusters
dat.long.std.ave.mode <- dat.long.std.ave %>%
  group_by(cluster, latitude) %>%
  summarise(mean_std_area = as.numeric(mean(std_area)),
            stdev_std_area = as.numeric(sd(std_area))) %>%
  ungroup()


g.mode <- ggplot(data = dat.long.std.ave.mode, aes(x = as.numeric(as.character(latitude)), y = mean_std_area))+
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = mean_std_area - stdev_std_area,
                  ymax = mean_std_area + stdev_std_area), alpha = 0.2) +
  facet_wrap(cluster ~ ., ncol = 3)

cluster.plots[[3]] <- g.mode

 
  plot.list.archtype <- list()
  
  for(cluster.number in 1:numberofclusters){
    
    mass.features <- dat.long.std.ave %>% filter(cluster == cluster.number) %>% 
      ungroup() %>% select(MassFeature_Column) %>% unique()
    num.mass.features <- length(mass.features$MassFeature_Column)
    
    knowns.dat <- dat.long.std %>%
      left_join(combo.wide %>% select(MassFeature_Column, Identification, Confidence)) %>%
      filter(Confidence ==  1) %>% 
      filter(MassFeature_Column %in% mass.features$MassFeature_Column) 
    num.knowns.mass.features <- length(unique((knowns.dat$Identification)))
    knowns.mass.features <- unique((knowns.dat$Identification))
    #three.rando.knowns <- sample(knowns.dat$MassFeature_Column, 3)
    
    culture.seen <- culture.dat.2.ave %>%
      filter(MassFeature_Column %in% mass.features$MassFeature_Column) %>%
      filter(!is.na(LogBioVol)) %>% 
      ungroup() %>% select(MassFeature_Column) %>% unique()
    
    # knowns <- dat.long.std %>%
    #   left_join(combo.wide %>% select(MassFeature_Column, Identification)) %>%
    #   filter(MassFeature_Column %in% three.rando.knowns) %>%
    #   mutate(latitude = as.factor(round(latitude, digits =0))) %>%
    #   group_by(latitude, MassFeature_Column, Identification) %>%
    #   summarize(std_area = mean(std_area, na.rm = TRUE),
    #             stdev_std_area = sd(std_area, na.rm = TRUE))
    
    #Plot of the mode
    mode.plot <- ggplot(data = dat.long.std.ave.mode %>% filter(cluster == cluster.number), 
                        aes(x = as.numeric(as.character(latitude)), y = mean_std_area))+
      geom_point() +
      geom_line() +
      geom_ribbon(aes(ymin = mean_std_area - stdev_std_area,
                      ymax = mean_std_area + stdev_std_area), alpha = 0.2) + 
      annotate("text", x = 30, y = 1.8, 
               label = c(paste("Cluster ", cluster.number, sep = "")),
               fontface =2, size =8) + 
      annotate("text", x = 30, y = 1.5, 
               label = c(paste(num.mass.features, " mass features", sep = "")),
               size =6) +
      labs(y= "Standardized peak area", x = "Latitude")
    
    #Pie chart of how many are IDd
    pie.dat.1 <- data.frame(
      group = c("Unknown", "Identified"),
      value = c(num.mass.features-num.knowns.mass.features, num.knowns.mass.features))
    pie.dat.2 <- data.frame(
      group = c("Not in culture collection", "In culture collection"),
      value = c(num.mass.features-length(culture.seen$MassFeature_Column), length(culture.seen$MassFeature_Column)))
    
    
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.title=element_text(size=14, face="bold")
      )
    
    pie.plot.1 <- ggplot(pie.dat.1, aes(x="", y=value, fill=group))+
      geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
      scale_fill_manual(values=c("#E69F00", "#999999")) +
      blank_theme +
      geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                    label = paste(value, "\n", group)), size=5)
    
    pie.plot.2 <- ggplot(pie.dat.2, aes(x="", y=value, fill=group))+
      geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
      scale_fill_manual(values=c("#E69F00", "#999999")) +
      blank_theme +
      geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                    label = paste(value, "\n", group)), size=5)
    
    pie.plots <- plot_grid(pie.plot.1, pie.plot.2, ncol = 1)
    
    #Put these three together
    mode.pies <- plot_grid(mode.plot, pie.plots, ncol = 2, rel_widths = c(1, 0.8), align = "h")
    
    plot.list.archtype[[cluster.number]] = mode.pies
  }
   
  cluster.plots[[4]] <- plot.list.archtype
  
  return(cluster.plots)
  
  }
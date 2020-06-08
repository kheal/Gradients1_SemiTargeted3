library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


#Get the plots from the dependent files
source("Figures/Code/ctd_plots.R")
source("Figures/Code/tile_plots_all.R")
source("Figures/Code/enviro_cluster_mode_plots.R")
source("Tables/Code/cluster_compounds.R")
source("Figures/Code/SampleMap.R")


#Make a multipaneled plot
g.MGL.all.sub <- plot_grid(g.MGL.ctd, g.tile.MGL, ncol = 1, align = "v", axis = "bt", 
                           rel_heights = c(1, 1), labels = c("A", "B"))
g.MGL.all.sub.2 <- plot_grid(g.MGL.all.sub, g.mode.MGL, ncol = 2, labels = c("", "C"))
g.MGL.all.sub.3 <- plot_grid(g.MGL.all.sub.2, g.cmps.MGL, ncol = 1,  rel_heights = c(1, 0.5), labels = c("", "D"))
g.MGL.all.sub.3
save_plot("Figures/Manuscript_figures/MGL_modeplots_multipanel.pdf", g.MGL.all.sub.3, base_height = 9, base_width = 8)

g.KM.all.sub <- plot_grid(g.KM.ctd, g.tile.KM, ncol = 1, align = "v", axis = "bt", rel_heights = c(1, 1), labels = c("A", "B"))
g.KM.all.sub.2 <- plot_grid(g.KM.all.sub, g.mode.KM, ncol = 2,  labels = c("", "C"))
g.KM.all.sub.3 <- plot_grid(g.KM.all.sub.2, g.cmps.KM, ncol = 1,  rel_heights = c(1, 0.5), labels = c("", "D"))
g.KM.all.sub.3
save_plot("Figures/Manuscript_figures/KM_modeplots_multipanel.pdf", g.KM.all.sub.3, base_height = 9, base_width = 8)

g.KOK.all.sub <- plot_grid(KOK.map, g.KOK.ctd, ncol = 1, rel_heights = c(1.8, 1), rel_widths = c(1, 1.8), labels = c("A", "C"))
g.KOK.all.sub.2 <- plot_grid(g.KOK.all.sub, g.mode.KOK, ncol = 2,  rel_widths = c(1, 1.3), labels = c("", "C"))
g.KOK.all.sub.2
g.KOK.all.sub.3 <- plot_grid(g.KOK.all.sub.2, g.cmps.KOK, ncol = 1,  rel_heights = c(1, 0.3), labels = c("", "D"))
g.KOK.all.sub.3
save_plot("Figures/Manuscript_figures/KOK_modeplots_multipanel.pdf", g.KOK.all.sub.3, base_height = 10, base_width = 8)


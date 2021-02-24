# To combine both NMDS plots into one for publication
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)

p1+ labs(subtitle:"A") + p2 +labs(subtitle = "B")


d3 <- (d + labs(tag = 'A'))/ (d.raw  +labs(tag = 'B'))
d3

save_plot("Figures/Manuscript_figures/NMDS_combo.pdf", d3, base_height = 8, base_width = 6)

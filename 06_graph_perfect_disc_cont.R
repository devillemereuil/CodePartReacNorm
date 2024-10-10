conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(tidyr::unpack)

# Libraries
library(tidyverse)
library(here)
library(patchwork)
library(grid)

# Some ggplot options
theme_set(theme_bw())
theme_update(text = element_text(family = "Linux Biolinum O", size = 20),
             plot.title = element_text(face = "bold", hjust = 0.5),
             axis.title.x = element_text(size = 24),
             axis.text.x = element_text(vjust = 0.5))


## ------- Loading the ggplot objects to plot together in the hybrid graph ----

p_part_disc <- readRDS(here("ggplot_part_perfect_disc_sub.rds"))
p_box_disc <- readRDS(here("ggplot_box_perfect_disc_sub.rds"))
p_part_cont <- readRDS(here("ggplot_part_perfect_cont_sub.rds"))
p_box_cont <- readRDS(here("ggplot_box_perfect_cont_sub.rds"))

ylab <- p_part_disc[["labels"]][["y"]]
p_part_disc[["labels"]][["y"]] <- p_part_cont[["labels"]][["y"]] <- ""

## --------------------------------------------------- Generating the plot ----

cairo_pdf(here("Figs/Sim_varpart_bias_hybrid.pdf"), height = 10, width = 15)
plot(
    (((p_part_disc + theme(axis.title.x = element_blank()) | p_box_disc) + 
       plot_layout(widths = c(4, 1))) / 
     ((p_part_cont | p_box_cont) + plot_layout(widths = c(4, 1)))) +
     plot_layout(heights = c(1.8, 1.2), axis_titles = "collect")
)
grid.draw(textGrob(ylab, x = 0.02, rot = 90,
                   gp = gpar(fontsize = 24, fontfamily = "Linux Biolinum O")))
dev.off()
